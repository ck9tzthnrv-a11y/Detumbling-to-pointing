% Esempio di chiamata valida (sostituire i valori con i dati reali)
q_pred    = [0.0234;-0.0422;0.0269;0.9985];        % quaternione predetto [4x1]
omega_pred= [-0.3595;-0.1770;-0.1183];          % velocita angolare predetta [3x1]
P_pred    = eye(6);           % covarianza predetta [6x6]
q_meas    = [0;0;0;1];        % quaternione misurato [4x1]
omega_meas= [0;0;0];          % velocita angolare misurata [3x1]
R         = eye(6);           % matrice di covarianza delle misure [6x6]

[q_post, omega_post, P_post] = mekf_correction_step(q_pred, omega_pred, P_pred, q_meas, omega_meas, R);
disp('q_post:'), disp(q_post')
disp('omega_post (deg/s):'), disp(omega_post'*180/pi)
disp('P_post diagonal:'), disp(diag(P_post)')

function [q_post, omega_post, P_post] = mekf_correction_step(q_pred, omega_pred, P_pred, q_meas, omega_meas, R)

% =========================================================================
% MEKF CORRECTION STEP - Aggiornamento stato e covarianza
% =========================================================================

% INPUT:
%   q_pred     : quaternione predetto [4x1], scalar-last, da mekf_prediction_step
%   omega_pred : velocita' angolare predetta [3x1] [rad/s], da mekf_prediction_step
%   P_pred     : covarianza predetta [6x6], da mekf_prediction_step
%   q_meas     : quaternione misurato dallo Star Tracker [4x1], scalar-last
%                (ASTRO APS: accuratezza < 1 arcsec across boresight, 1-sigma)
%   omega_meas : velocita' angolare misurata dai giroscopi [3x1] [rad/s]
%                (Astrix 1090: ARW < 0.005 deg/sqrt(h), bias < 0.01 deg/h)
%   R          : matrice di covarianza delle misure [6x6]
%                Struttura:
%                  std_phidq    = 1 deg (Star Tracker, 1-sigma)
%                  Rattitude    = (std_phidq*pi/180)^2 / 12 * eye(3)
%                  Rangvel      = Sigma_omegameas^2 * eye(3)
%                  R = blkdiag(Rattitude, Rangvel)
%
% OUTPUT:
%   q_post     : quaternione corretto e rinormalizzato [4x1], scalar-last
%   omega_post : velocita' angolare corretta [3x1] [rad/s]
%   P_post     : covarianza aggiornata [6x6], simmetrica e definita positiva

% =========================================================================
% 1. INNOVAZIONE (RESIDUO)
% =========================================================================

% --- Innovazione di assetto ---
% Calcola il quaternione errore tra misura e predizione:
%   dq = q_meas * q_pred^{-1}
% Per un quaternione unitario: q^{-1} = q* = [-qv; qs]
q_pred_conj = [-q_pred(1:3); q_pred(4)];            % coniugato = inverso
dq_inn = q_product(q_meas, q_pred_conj);             % quaternione errore

% Correzione emisferica: impone dq_s >= 0 per evitare salti di segno
% nell'innovazione che causerebbero instabilita' durante la convergenza.
if dq_inn(4) < 0
    dq_inn = -dq_inn;
end

% Innovazione in termini di vettore di Gibbs: dg = dq_v / dq_s
% (approssimazione: dg ≈ dθ/2 per piccoli angoli)
% La distinzione e' importante durante la convergenza iniziale.
delta_g_inn = dq_inn(1:3) ./ dq_inn(4);             % vettore di Gibbs [3x1]

% --- Innovazione di velocita' angolare ---
delta_omega_inn = omega_meas - omega_pred;           % residuo omega [3x1]

% --- Vettore di innovazione completo ---
innovation = [delta_g_inn; delta_omega_inn];         % [6x1]

% =========================================================================
% 2. GUADAGNO DI KALMAN
% =========================================================================
% H = I(6x6): il modello di misura e' diretto (stato = misura).
H = eye(6);

% Matrice di innovazione S = H*P*H' + R = P_pred + R (poiche' H=I)
S = P_pred + R;                                      % [6x6]

% Guadagno ottimale: K = P * H' * S^{-1} = P_pred * S^{-1}
K = P_pred * (S \ eye(6));                           % [6x6], equivale a P/S

% =========================================================================
% 3. AGGIORNAMENTO DELLO STATO D'ERRORE
% =========================================================================
% Correzione nello spazio degli errori: delta_x = K * innovation
delta_x = K * innovation;                            % [6x1]

dg_up  = delta_x(1:3);   % errore d'assetto stimato (Gibbs vector) [3x1]
dw_up  = delta_x(4:6);   % errore di velocita' angolare stimato [3x1]

% =========================================================================
% 4. RESET (iniezione dell'errore nello stato globale)
% =========================================================================

% --- Aggiornamento moltiplicativo del quaternione ---
% Converte il vettore di Gibbs dg in quaternione d'errore normalizzato.
% Formula: dq_k = [dg; 1] / sqrt(1 + ||dg||^2)
dq_up = [dg_up; 1] / sqrt(1 + dg_up' * dg_up);     % quaternione d'errore [4x1]

% Composizione moltiplicativa: q_post = dq_up * q_pred
q_post = q_product(dq_up, q_pred);
q_post = q_post / norm(q_post);                      % rinormalizzazione di sicurezza

% --- Aggiornamento additivo della velocita' angolare ---
omega_post = omega_pred + dw_up;                     % [3x1], rad/s

% =========================================================================
% 5. AGGIORNAMENTO DELLA COVARIANZA
% =========================================================================
% Forma di Joseph: P = (I-KH)*P*(I-KH)' + K*R*K'
% Garantisce la definita positività di P anche con errori numerici,
% necessario per #codegen su hardware embedded (GR740 OBC).
% La forma standard (I-KH)*P e' numericamente instabile se K non e' ottimale.
IKH   = eye(6) - K * H;                             % [6x6]
P_post = IKH * P_pred * IKH' + K * R * K';          % forma di Joseph

% Simmetrizzazione finale per robustezza numerica
P_post = (P_post + P_post') / 2;

end
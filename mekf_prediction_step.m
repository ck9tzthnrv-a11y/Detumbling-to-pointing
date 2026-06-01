function [q_pred, omega_pred, P_pred] = mekf_prediction_step(q_prev, omega_prev, P_prev, dt, J, Qangvelfactor)
% MEKF PREDICTION STEP - Propagazione stato e covarianza

% INPUT:
%   q_prev        : quaternione [4x1] al passo k-1, scalar-last [q1 q2 q3 q4]
%                   All'inizializzazione: q_prev = q_opt da solve_wahba_q_method
%   omega_prev    : velocita' angolare [3x1] al passo k-1 [rad/s]
%   P_prev        : covarianza [6x6] al passo k-1
%                   All'inizializzazione: assemblare come
%                     P_prev = [P_wahba, zeros(3,3); zeros(3,3), sigma_omega^2*eye(3)]
%                   con P_wahba da solve_wahba_q_method e
%                   sigma_omega = 2e-4 rad/s (Astrix 1090, datasheet)
%   dt            : passo di integrazione [s]
%   J             : matrice d'inerzia 3x3 [kg*m^2]
%                   Per il progetto: diag([7600, 8700, 5100])
%   Qangvelfactor : PSD del rumore angolare [rad^2/s^3], scalare
%                   Stima da Astrix 1090 ARW = 0.005 deg/sqrt(h):
%                   ARW_rad = 0.005/sqrt(3600)*pi/180 = 1.45e-6 rad/sqrt(s)
%                   Qangvelfactor = ARW_rad^2 = ~2.1e-12 rad^2/s^3
%                   In pratica nel progetto si usa 1e-9 (margine conservativo)
%
% OUTPUT:
%   q_pred     : quaternione propagato [4x1], scalar-last, rinormalizzato
%   omega_pred : velocita' angolare propagata [3x1] [rad/s]
%   P_pred     : covarianza propagata [6x6]

%% 1. PROPAGAZIONE QUATERNIONE
% Usa propagate_quaternion (funzione del progetto) che implementa
% l'integrazione della cinematica q_dot = (1/2)*Omega(omega)*q
q_pred = propagate_quaternion(q_prev, omega_prev, dt);
% propagate_quaternion include gia' la rinormalizzazione interna.

%% 2. PROPAGAZIONE VELOCITA' ANGOLARE (Equazioni di Eulero senza coppia esterna)
% omega_dot = J^{-1} * (-(omega x J*omega))
% La coppia esterna e' gestita nel correction step o come disturbance
% nel Q_proc, non nella predizione nominale.
skew_omega = [0,            -omega_prev(3),  omega_prev(2);
    omega_prev(3),  0,            -omega_prev(1);
    -omega_prev(2),  omega_prev(1),  0           ];

omega_dot  = inv(J) * (-(skew_omega * J * omega_prev));
omega_pred = omega_prev + omega_dot .* dt;

%% 3. JACOBIANO DEL PROCESSO F (6x6)
% Linearizzazione attorno allo stato corrente, coerente con JacobProcess
Jomega     = J * omega_prev;
skew_Jomega = [0,          -Jomega(3),  Jomega(2);
    Jomega(3),   0,         -Jomega(1);
    -Jomega(2),   Jomega(1),  0         ];

F = [-skew_omega,    0.5*eye(3);
    zeros(3,3),    inv(J)*(-skew_omega*J + skew_Jomega)];

%% 4. MATRICE DI TRANSIZIONE Phi (esponenziale di matrice)
% expm e' piu' accurata di eye(6)+F*dt per dt non trascurabili.
Phi = expm(F .* dt);

%% 5. Q DISCRETA (Metodo di Van Loan)
% G: matrice di ingresso del rumore - il rumore di processo agisce SOLO
% sulle 3 componenti di velocita' angolare (non sull'attitude direttamente).
G = [zeros(3,3); eye(3)];

% Van Loan: trasforma l'integrale della covarianza discreta in un esponenziale
% Qd = int_0^dt [Phi(tau)*G*q*G'*Phi(tau)'] dtau
M = [-F,                  G * Qangvelfactor * G';
    zeros(6),           F'                   ] .* dt;

E   = expm(M);
Phi_vl = E(7:12, 7:12)';     % recupera Phi dal blocco in basso a destra
Qd  = Phi_vl * E(1:6, 7:12); % covarianza discreta

%% 6. PROPAGAZIONE COVARIANZA
P_pred = Phi * P_prev * Phi' + Qd;

% Simmetrizzazione numerica (previene deriva per accumulo di errori float)
P_pred = (P_pred + P_pred') / 2;

end
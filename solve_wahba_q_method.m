function [q_opt, P_wahba] = solve_wahba_q_method(B_matrix, R_matrix, sigma_vec)
% Risolutore Wahba Q-method (scalar-last quaternion)

% INPUT:
%   B_matrix : 3xN vettori misurati in body frame (colonne, non normalizzati)
%   R_matrix : 3xN vettori di riferimento in reference frame (colonne)
%   sigma_vec: 1xN deviazioni standard angolari delle misure [gradi],
%   ovvero quanto sono precisi i sensori

% OUTPUT:
%   q_opt   : quaternione ottimale [q1 q2 q3 q4], scalar-last
%   P_wahba : covarianza 3x3 dell'errore in dtheta [rad^2]

if nargin == 0
    % --- BLOCCO DI AUTODIAGNOSTICA / UNIT TESTING ---
    % Best practice: permette di runnare lo script a vuoto per verificare che funzioni.
    % Test di default coerente col progetto ADCS (2 sensori, es. Sun Sensor e Magnetometro)
    B_matrix = [0; -1; 0];           % es. sun sensor body
    B_matrix(:,2) = [0.6; 0; 0.8];   % es. magnetometro body
    R_matrix = [0; -1; 0];
    R_matrix(:,2) = [0.6; 0; 0.8];
    sigma_vec = [0.01, 0.5];         % gradi: sun sensor, magnetometro
    [q_opt, P_wahba] = solve_wahba_q_method(B_matrix, R_matrix, sigma_vec);
    disp('q_opt:'); disp(q_opt');
    disp('P_wahba sqrt diag (deg):');
    disp(sqrt(diag(P_wahba))' * 180/pi);
    return
end

N = size(B_matrix, 2);
sigma_rad = deg2rad(sigma_vec);   % converti in radianti

% --- CALCOLO PESI NORMALIZZATI ---
% L'algoritmo deve fidarsi di più dei sensori precisi.
% sigmatot rappresenta la deviazione standard equivalente globale del set.
sigmatot = 1 / sqrt(sum(1 ./ sigma_rad.^2));
a = sigmatot^2 ./ sigma_rad.^2;   % Pesi normalizzati: sum(a) = 1

% --- Attitude Profile Matrix B ---
Bmat = zeros(3,3);
z    = zeros(3,1);
for i = 1:N
    % Normalizzazione di sicurezza per evitare errori da sensori/modelli
    b_i = B_matrix(:,i) / norm(B_matrix(:,i));
    r_i = R_matrix(:,i) / norm(R_matrix(:,i));
    % Bmat: Matrice 3x3 che condensa tutte le info direzionali pesate.
    Bmat = Bmat + a(i) * (b_i * r_i');
    % z: Vettore che cattura l'errore fisico di disallineamento incrociato.
    z    = z    + a(i) * cross(b_i, r_i);
end

% --- MATRICE DI DAVENPORT (K) ---
% Davenport trasforma il problema di minimizzazione di Wahba in un
% problema agli autovalori sui quaternioni.
S     = Bmat + Bmat';
sigma = trace(Bmat);
K = [ S - sigma*eye(3),  z ;
    z'              ,  sigma ]; % Struttura valida per convenzione SCALAR-LAST

% --- SOLUZIONE AGLI AUTOVALORI (IL Q-METHOD) ---
[V, D] = eig(K);
% Il quaternione che minimizza la loss function è l'autovettore
% corrispondente all'autovalore MASSIMO della matrice K.
[~, idx] = max(diag(D));
q_opt = V(:, idx);

% --- SANITY CHECK: AMBIGUITÀ EMISFERICA ---
% I quaternioni q e -q descrivono la stessa rotazione fisica.
% Forziamo la parte scalare positiva per evitare fastidiosi "salti di segno" 
% (unwinding) che manderebbero in instabilità la propagazione nel MEKF.
if q_opt(4) < 0
    q_opt = -q_opt;
end

% --- COVARIANZA IN dtheta [rad^2] (FISHER INFORMATION MATRIX) ---
% Generiamo il dato sul rumore della misura d'assetto da passare al Kalman Filter.
sum_weighted_outer = zeros(3,3);
for i = 1:N
    b_i = B_matrix(:,i) / norm(B_matrix(:,i));
    sum_weighted_outer = sum_weighted_outer + a(i) * (b_i * b_i');
end

% Info_matrix quantifica la bontà geometrica delle misurazioni attuali.
Info_matrix = eye(3) - sum_weighted_outer;

% --- CONTROLLO DI ROBUSTEZZA ---
% Se i vettori misurati sono paralleli (es. Sole e Campo Magnetico allineati), 
% si perde un grado di libertà e la matrice diventa quasi-singolare.
if rcond(Info_matrix) < 1e-12
    warning('Info_matrix quasi-singolare: configurazione geometrica degenere.');
    % Evita il crash di Simulink calcolando la pseudo-inversa (pinv)
    P_wahba = pinv(Info_matrix) * sigmatot^2;
else
    % Condizione operativa nominale
    P_wahba = sigmatot^2 * inv(Info_matrix);
end
end
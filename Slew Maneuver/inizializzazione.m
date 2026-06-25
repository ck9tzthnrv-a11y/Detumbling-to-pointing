%% SCRIPT DI INIZIALIZZAZIONE ADCS - SATELLITE 1 TONNELLATA
clear; clc;

%% 1. PARAMETRI DEL CORPO RIGIDO E HARDWARE
% Matrice di Inerzia (Worst-case scenario J = 5000)
J = [3190 0 0; 0 3790 0; 0 0 5672]; % [kg*m^2] Inerzie sui 3 assi principali

% Limiti Hardware delle Ruote di Reazione (Collins RSI 12-75/60)
tau_m = 0.106;          % Coppia massima erogabile per ruota [Nm]
h_max = 12;             % Capacità di momento per ruota [Nms]
theta_deg = 45;
% 2. PARAMETRI DEL CONTROLLORE DI KJELLBERG
J_max = max(diag(J));
theta_dot_m = (2 * h_max) / J_max; % Rateo massimo [rad/s]

% --- NUOVO TUNING BASATO SUI LIMITI HARDWARE ---
% Imponiamo che la richiesta massima di coppia K_max non superi tau_m
% K_max = c * theta_dot_m * J_max  --->  c = tau_m / (theta_dot_m * J_max)

% c = tau_m / (theta_dot_m * J_max);   % "Freno" derivativo tarato sull'hardware
% kappa = (c^2)/2 ;                   % "Molla" prop. per garantire smorzamento critico (zeta = 0.707)

% Tuning del controllore PD per Inseguimento Orbitale (Tracking)
omega_n = 0.05;         % Frequenza naturale desiderata [rad/s] (es. 0.05)
zeta = 1.0;             % Smorzamento Critico per tracking senza oscillazioni

kappa = omega_n^2;      % Guadagno "molla" (corrisponde a Kp diviso l'inerzia)
c = 2 * zeta * omega_n; % Guadagno "freno" (corrisponde a Kd diviso l'inerzia)


%% 4. PARAMETRI ORBITALI (Per calcolo disturbi ambientali)
% Dati per una tipica orbita bassa (LEO) - Esempio tratto dai tuoi file
R_E = 6378.1363e3;      % Raggio della Terra [m]
mu = 3.986e14;          % Parametro gravitazionale terrestre [m^3/s^2]
alt = 600e3;            % Altitudine [m]

a = R_E + alt;          % Semiasse maggiore [m]
ecc = 1e-10;            % Eccentricità (Quasi circolare)
inc = 15.3 * (pi/180);    % Inclinazione  [rad]
RAAN = 0;               % Right Ascension of Ascending Node [rad]
AOP = 0;                % Argomento del perigeo  [rad]
ni = 0;                 % Anomalia vera iniziale [rad]

omega_orb = sqrt(mu / a^3); % Velocità angolare orbitale (Mean motion) [rad/s]
%% 3. STATO INIZIALE E ASSETTO BERSAGLIO (ALLINEAMENTO CON ORF/LVLH)
% Stato Iniziale Vero del satellite (Fermo e allineato all'inerziale ECI)
q_true_0 = [ -0.7125
   -0.2789
    0.2955
    0.5720];
w_true_0 = [-0.0191;-0.0438;-0.237]*pi/180;
state_0 = [q_true_0; w_true_0]; % Vettore 7x1 per l'Integratore

% Calcolo Posizione (r) e Velocità (v) in ECI all'istante iniziale t=0
% Assumendo satellite all'equatore (Anomalia = 0, RAAN = 0, AOP = 0)
v_orb = sqrt(mu / a); 
r_ECI = [a; 0; 0];
v_ECI = [0; v_orb * cos(inc); v_orb * sin(inc)];

% COSTRUZIONE DEI VERSORI DELL'ORBITAL REFERENCE FRAME (LVLH)
z_o = -r_ECI / norm(r_ECI);                 % Asse Z: Nadir (verso la Terra)
h_orb = cross(r_ECI, v_ECI);
y_o = -h_orb / norm(h_orb);                 % Asse Y: Normale all'orbita (negativa)
x_o = cross(y_o, z_o);                      % Asse X: Direzione del moto

% Matrice di rotazione da ECI a ORF
A_IO = [x_o, y_o, z_o]; % Trasforma da ORF a ECI
A_OI = A_IO';           % Trasforma da ECI a ORF
% Estrazione ROBUSTA del quaternione q_target dalla matrice A_OI (DCM)
% Basato sul metodo di Markley/Shepperd per evitare divisioni per zero
tr = trace(A_OI);

% Troviamo il valore massimo tra la traccia e gli elementi sulla diagonale
[~, max_idx] = max([tr, A_OI(1,1), A_OI(2,2), A_OI(3,3)]);

switch max_idx
    case 1 % La traccia è dominante (rotazioni piccole/medie)
        q4 = 0.5 * sqrt(1 + tr);
        q1 = (A_OI(2,3) - A_OI(3,2)) / (4*q4);
        q2 = (A_OI(3,1) - A_OI(1,3)) / (4*q4);
        q3 = (A_OI(1,2) - A_OI(2,1)) / (4*q4);
        
    case 2 % L'elemento A_OI(1,1) è dominante (rotazione ampia attorno X)
        q1 = 0.5 * sqrt(1 + A_OI(1,1) - A_OI(2,2) - A_OI(3,3));
        q4 = (A_OI(2,3) - A_OI(3,2)) / (4*q1);
        q2 = (A_OI(1,2) + A_OI(2,1)) / (4*q1);
        q3 = (A_OI(3,1) + A_OI(1,3)) / (4*q1);
        
    case 3 % L'elemento A_OI(2,2) è dominante (rotazione ampia attorno Y)
        q2 = 0.5 * sqrt(1 + A_OI(2,2) - A_OI(1,1) - A_OI(3,3));
        q4 = (A_OI(3,1) - A_OI(1,3)) / (4*q2);
        q1 = (A_OI(1,2) + A_OI(2,1)) / (4*q2);
        q3 = (A_OI(2,3) + A_OI(3,2)) / (4*q2);
        
    case 4 % L'elemento A_OI(3,3) è dominante (rotazione ampia attorno Z)
        q3 = 0.5 * sqrt(1 + A_OI(3,3) - A_OI(1,1) - A_OI(2,2));
        q4 = (A_OI(1,2) - A_OI(2,1)) / (4*q3);
        q1 = (A_OI(3,1) + A_OI(1,3)) / (4*q3);
        q2 = (A_OI(2,3) + A_OI(3,2)) / (4*q3);
end
% Quaternione target in convenzione Scalar-Last
q_target = [q1; q2; q3; q4]; 

disp('Inizializzazione completata! Puoi avviare Simulink.');
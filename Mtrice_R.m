function R = mekf_r_matrix_block()
%% =========================================================
%% COSTRUZIONE DELLA MATRICE R PER IL MEKF CORRECTION STEP
%% =========================================================
% Coerente con i datasheet dei sensori del progetto.

%% --- Blocco attitude (Star Tracker ASTRO APS) ---
% Fonte: Data_Sheet_ASTRO_APS_StarSensor.pdf
% Accuratezza across boresight: < 1 arcsec [1-sigma]
% Accuratezza boresight:        < 8 arcsec [1-sigma]
% Nel file di riferimento si usa std_phidq = 1 deg come valore
% conservativo (include margini di sistema e misallineamento).

std_phidq   = 1;                           % [deg], 1-sigma, valore di progetto
std_phidq_r = std_phidq * pi/180;          % [rad]

% Rattitudefactor: varianza su ciascuna componente di dg
% Derivazione: Var(dg_i) = sigma_phi^2 / 4 (da dg≈dθ/2)
%              diviso 3 per distribuzione isotropica su sfera
%              --> fattore totale = 1/12
Rattitudefactor = std_phidq_r^2 / 12;     % [rad^2], scalare

R_att = eye(3) * Rattitudefactor;          % [3x3]

%% --- Blocco velocità angolare (Astrix 1090) ---
% Fonte: satbase-astrix-1090-IMU.pdf
% ARW = 0.005 deg/sqrt(h) [1-sigma]
% Bias stability (1h) = 0.01 deg/h
% Nel file di riferimento: Sigma_omegameas = 2e-4 rad/s
% Conversione ARW in sigma istantanea a frequenza f_meas [Hz]:
%   sigma_omega = ARW_rad / sqrt(1/f_meas)
%               = ARW_rad * sqrt(f_meas)
% Con f_meas = 2 Hz (come nel file di riferimento, measfrequency=2):
%   ARW_rad = 0.005/sqrt(3600) * pi/180 = 1.454e-6 rad/sqrt(s)
%   sigma_omega = 1.454e-6 * sqrt(2) = 2.06e-6 rad/s
% Il file di riferimento usa il valore conservativo 2e-4 rad/s
% (include deriva e altri effetti a lungo termine).

Sigma_omegameas = 2e-4;                    % [rad/s], da file di riferimento
Rangvelfactor   = Sigma_omegameas^2;       % [rad^2/s^2]

R_omega = eye(3) * Rangvelfactor;          % [3x3]

%% --- Assemblaggio R completa ---
R = [R_att,         zeros(3,3);
    zeros(3,3),    R_omega];              % [6x6]
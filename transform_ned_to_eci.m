function mag_ref_eci = transform_ned_to_eci(mag_ned, lat_deg, lon_deg, gmst_deg)
% TRASFORMAZIONE CAMPO MAGNETICO: NED -> ECEF -> ECI
% INPUTS:
%   mag_ned  : Collegare DIRETTAMENTE al PRIMO OUTPUT del blocco IGRF 13
%   lat_deg  : Latitudine geodetica istantanea [gradi] (dalla propagazione orbitale)
%   lon_deg  : Longitudine geodetica istantanea [gradi] (dalla propagazione orbitale)
%   gmst_deg : Greenwich Mean Sidereal Time [gradi] (da Julian Date to GMST block)
% OUTPUT:
%   mag_ref_eci: Vettore direzione 3x1 normalizzato in ECI (per Wahba R_matrix)

% --- SICUREZZA INGRESSO (Robusta per Simulink) ---
% Forza il vettore in ingresso a essere una colonna 3x1,
% indipendentemente se Simulink lo invia come riga (1x3) o colonna (3x1).
mag_ned = mag_ned(:); 

% Conversione in radianti per le funzioni trigonometriche
lat = deg2rad(lat_deg);
lon = deg2rad(lon_deg);
gmst = deg2rad(gmst_deg);

% --- 1. MATRICE DI ROTAZIONE: NED -> ECEF ---
% Basata sulle coordinate geodetiche del satellite
R_ned_to_ecef = [ -sin(lat)*cos(lon), -sin(lon), -cos(lat)*cos(lon);
    -sin(lat)*sin(lon),  cos(lon), -cos(lat)*sin(lon);
    cos(lat),           0,        -sin(lat) ];

% --- 2. MATRICE DI ROTAZIONE: ECEF -> ECI ---
% Basata sulla rotazione terrestre istantanea (GMST)
R_ecef_to_eci = [ cos(gmst), -sin(gmst), 0;
    sin(gmst),  cos(gmst), 0;
    0,          0,         1 ];

% --- 3. TRASFORMAZIONE COMPLETA ---
% v_eci = R(ECEF->ECI) * R(NED->ECEF) * v_ned
mag_eci_unnorm = R_ecef_to_eci * R_ned_to_ecef * mag_ned;

% --- 4. NORMALIZZAZIONE (Essenziale per Wahba) ---
% Eliminiamo l'intensità (nT) e teniamo solo la direzione
norm_val = norm(mag_eci_unnorm);
if norm_val < 1e-9 % Protezione contro divisione per zero
    mag_ref_eci = [0; 0; 0]; 
else
    mag_ref_eci = mag_eci_unnorm / norm_val;
end
end
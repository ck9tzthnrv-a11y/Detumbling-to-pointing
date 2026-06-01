function P0 = matrice_P_0(P_wahba)
%#codegen
% --- Inizializzazione Matrice di Covarianza P0 per MEKF ---
% INPUT:
%   P_wahba : Matrice di covarianza 3x3 dell'assetto calcolata da Wahba al tempo t=0
% OUTPUT:
%   P0      : Matrice di covarianza globale 6x6 iniziale per il filtro

    % --- SICUREZZA INGRESSO ---
    % Garantiamo che Simulink non dia errori se la matrice arriva "piatta"
    if ~isequal(size(P_wahba), [3, 3])
        P_wahba = reshape(P_wahba, 3, 3);
    end
  
    % 1. Incertezza Velocità Angolare (Giroscopi)
    % Valore da datasheet Astrix 1090 (2e-4 rad/s) 
    sigma_omega_init = 2e-4;
    P_omega_init = (sigma_omega_init^2) * eye(3);

    % 2. Assemblaggio Matrice P0 (6x6)
    % Struttura a blocchi diagonali:
    % In alto a sinistra: incertezza sull'assetto (P_wahba)
    % In basso a destra: incertezza sulla velocità angolare (P_omega_init)
    P0 = [P_wahba,       zeros(3,3);
          zeros(3,3),    P_omega_init]
end
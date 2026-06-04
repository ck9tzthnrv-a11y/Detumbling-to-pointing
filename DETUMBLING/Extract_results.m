%% Estrazione condizioni iniziali da Simulink detumbling
q0_est = out.q_true(end-3:end);
w0_est = out.w_true(end-2:end);
init_state_est = [q0_est; w0_est];   % [7×1]



disp('Quaternione finale:')
disp(q0_est)
disp(['Norma q: ', num2str(norm(q0_est))])
disp('Velocità angolare finale [deg/s]:')
disp(w0_est * 180/pi)

% Salva TUTTE le variabili necessarie
save('detumbling_results.mat', ...
    'q0_est', 'w0_est', 'init_state_est', ...
    'noise_power_ARW', 'noise_power_RRW', 'Ts_gyro', ...
    'sample_f', 'ARW', 'bias_instab', 'PSD_ARW', 'PSD_bias','Sat_Inertia', 'J');
disp('Condizioni iniziali salvate in detumbling_results.mat')
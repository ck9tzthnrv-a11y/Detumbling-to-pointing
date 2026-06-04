clear all;
close all;
clc;

%% Simulation Settings
t_sim = 30000; 

%% Initializer
init_attitude_eul = [10.1, 0.1, 0.1];
init_attitude_rad = convang(init_attitude_eul,'deg','rad');
init_angvel = [3 3 3];
init_angvel_rad = convangvel(init_angvel,'deg/s','rad/s');
init_quat = eul2quat(init_attitude_rad);
init_quat = [init_quat(2:4), init_quat(1)];
init_state = [init_quat'; init_angvel_rad'];

desiredattitude = [0 0 0];
desiredquaternion = eul2quat(desiredattitude);
desiredquaternionSL = [desiredquaternion(2:4), desiredquaternion(1)];
attitudeprofile = timeseries(repmat(desiredquaternionSL, t_sim, 1));

%% Satellite Characteristics
Sat_Inertia = diag([7600, 8700, 5100]);
J = Sat_Inertia;

%% k_bdot
I_min = 5100;
alt = 600; 
mu_earth = 398600; 
R_earth = 6378.1363; 
a = R_earth + alt;
omega_0 = sqrt(mu_earth / a^3);
inc_m = 15.3 * (pi/180);
B_mean = 30000 * 1e-9;
k_bdot = (2 * omega_0 * (1 + sin(inc_m)) * I_min) / (B_mean^2);
disp(['k_bdot: ', num2str(k_bdot,'%e')]);

%% Gyroscope Parameters
sample_f = 10;
Ts_gyro  = 1/sample_f;
ARW      = 0.015;
PSD_ARW  = ARW * (pi/180) / 60;
noise_power_ARW = PSD_ARW^2;
bias_instab = 0.08;
PSD_bias    = bias_instab * (pi/180) / 3600;
noise_power_RRW = PSD_bias^2;


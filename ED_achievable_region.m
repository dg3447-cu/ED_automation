% -------------------------------------------------------------------------

% This script computes the closed capacity region for an envelope detector.

% The idea here is to bound all of the possible ED output voltage vs. # stages (N)
% into a closed region. This is useful for design automation of EDs. The
% idea is that if we write an optimization algorithm to design the optimal
% ED given a set of input parameters, the search space for the optimum ED
% will be closed, making it possible to find a global optimum ED. You can
% also use this script to manually choose how many stages to make your ED
% based off of sensitivity, BER, and data rate requirements.

% -------------------------------------------------------------------------

% Set simulation parameters

% RF parameters
P_RF = -60; % Target RX sensitivity (dBm)
RS = 3; % Source / antenna resistance (Ohms)
BER = 1e-3; % Desired BER
BW_BB = 1e3; % Baseband signal BW in Hz = Data rate

% Passive voltage gain from matching
Av = 30; % In dB

% Variance from simulation results after fabricationg the detector
fab_variation = 5; % Percent

% Detector Parameters
RD = 1e6; % Diode resistance (Ohms)
CC = 70e-15; % Coupling capacitance (F)

% Device parameters
n = 1.5; % Subthreshold slope constant (between 1 and 2)
Vt = 26e-3; % Thermal voltage
k = 1.38e-23; % Boltzmann's constant
T = 298; % Room temperature (K)

% -------------------------------------------------------------------------

% Set all of the bounds on the capacity region for ED output voltage vs. # of stages (N)

% Conversions
Av = 10 ^ (Av / 10);
P_RF_raw = P_RF;
P_RF = 10 ^ (P_RF / 10) * (10 ^ -3);

% Add a buffer for PVT variation
fab_variation = fab_variation / 100;

% Compute the upper bound for the number of allowed stages in the ED
r = 8 / (2.2 * RD * CC * BW_BB);
discriminant = sqrt(1 + r);
N_UPPER_BOUND = 0.25 * discriminant;
N_UPPER_BOUND = N_UPPER_BOUND * (1 - fab_variation);

% Compute the lower bound for the number of allowed stages in the ED
p_0 = exp(-1 / (pi + 2)) * sqrt(pi / (pi + 2));
num_lower_coeff = 5 * k * T;
Q_1_lower = -2 + 4 * ((0.25 + 0.5 * (BER / p_0)) ^ 0.5);
num_ln = log(Q_1_lower * p_0);
num_ln = num_ln * num_ln;
num_lower = num_lower_coeff * num_ln * (Vt ^ 2);
den_lower = (((1 / n) - 0.5) ^ 2) * (Av * P_RF * RS) ^ 2 * CC;
N_LOWER_BOUND = num_lower / den_lower;
N_LOWER_BOUND = N_LOWER_BOUND ^ (1 / 3);
N_LOWER_BOUND = N_LOWER_BOUND * (1 + fab_variation);

% Compute upper and lower bounds for ED voltage as a function of N
n_limits = 1:1:2; % Subthreshold slope constant is between 1 and 2
k_limits = (1 / (2 * Vt)) * ((1 ./ (n_limits)) - 0.5); % Min and Max conversion gain
N_values = 0:1:(floor(N_UPPER_BOUND * 1.2));
V_out_ED_MIN = (Av .^ 2) .* k_limits(2) .* N_values .* (P_RF * RS);
V_out_ED_MAX = (Av .^ 2) .* k_limits(1) .* N_values .* (P_RF * RS);

% Sanity check: if N_LOWER_BOUND > N_UPPER_BOUND, this indicates that the
% given RF parameters make it impossible to design an ED that meets the
% desired BER requirement.
if N_LOWER_BOUND > N_UPPER_BOUND
    fprintf("No ED exists for the presented BER requirement given the current simulation conditions.\n")
end

% Plot the capacity region
if N_LOWER_BOUND < N_UPPER_BOUND
    y_limits = [min(V_out_ED_MIN), max(V_out_ED_MAX)];
    figure; hold on
    
    % --- Highlight the parallelogram first (darker + more opaque)
    x_par = [N_LOWER_BOUND, N_UPPER_BOUND, N_UPPER_BOUND, N_LOWER_BOUND];
    y_par = [ ...
        (Av^2) * k_limits(2) * N_LOWER_BOUND * (P_RF * RS), ... % bottom-left
        (Av^2) * k_limits(2) * N_UPPER_BOUND * (P_RF * RS), ... % bottom-right
        (Av^2) * k_limits(1) * N_UPPER_BOUND * (P_RF * RS), ... % top-right
        (Av^2) * k_limits(1) * N_LOWER_BOUND * (P_RF * RS) ...  % top-left
    ];
    fill(x_par, y_par, [0.3 0.4 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none')
    
    % --- Then draw boundary lines on top
    plot(N_values, V_out_ED_MIN, 'b', 'LineWidth', 1.5)
    plot(N_values, V_out_ED_MAX, 'b', 'LineWidth', 1.5)
    plot([N_LOWER_BOUND N_LOWER_BOUND], y_limits, 'r', 'LineWidth', 1.5)
    plot([N_UPPER_BOUND N_UPPER_BOUND], y_limits, 'r', 'LineWidth', 1.5)
    
    xlabel("Number of detector stages (N)")
    ylabel("ED Output Voltage (V)")
    title("Closed Capacity Region for Envelope Detector")
    legend("Capacity Region", "V_{out, min}", "V_{out, max}", "N_{lower bound}", "N_{upper bound}", 'Location', 'northwest')
    grid on
    hold off
end

% -------------------------------------------------------------------------

% Plot three three system requirements again each other:
% Sensitivity = minimum RF input power required to achieve the given input parameters
% Data Rate
% BER

% Remember that these equations are derived off of bounds, so the provided
% plots are bounding regions where the curve itself gives the best case scenario.
% The third parameter that is not being swept is fixed at whatever the
% main simulation input is set to.

res = 1000;
figure;

% --- Plot bounding data rate vs. BER ---
BER_values = logspace(log10(0.00001), log10(0.1), res);

Q_1_lower_sweep = -2 + 4 * ((0.25 + 0.5 * (BER_values / p_0)) .^ 0.5);
num_ln_sweep = log(Q_1_lower_sweep * p_0);
num_ln_sweep = num_ln_sweep .* num_ln_sweep;
num_lower_sweep = num_lower_coeff * num_ln_sweep * (Vt ^ 2);

den_lower_sweep = (((1 / n) - 0.5) ^ 2) * (Av * P_RF * RS) ^ 2 * CC;
phi = num_lower_sweep / den_lower_sweep;
BER_function = 16 * ((phi) .^ (2 / 3)) - 1;
DataRate_values = 8 ./ (2.2 * RD * CC * BER_function);

subplot(2, 2, 1)
semilogx(BER_values, DataRate_values / (10 ^ 3)); grid;
xlabel("BER"); ylabel("Data Rate (kbps)"); title("Data Rate vs. BER Boundary");
legend(sprintf('Sensitivity (dBm) = %i', P_RF_raw))

% --- Plot bounding sensitivity vs. BER ---
num_coeff = 8000 * Vt * sqrt(5 * k * T);
num_ln = log(p_0 * (-2 + 4 * sqrt(0.25 + (BER_values / (2 * p_0)))));
num = num_coeff * num_ln;

den_coeff = Av * RS * sqrt(CC) * ((1 / n) - 0.5);
den_up_contr = 1 + (8 ./ (2.2 * RD * CC * BW_BB));
den_up_contr = (den_up_contr) .^ (3 / 4);
den = den_coeff .* den_up_contr;

P_RF_values_BER = abs(num ./ den);
P_RF_values_BER = 10 * log10(P_RF_values_BER);

subplot(2, 2, 2)
semilogx(BER_values, P_RF_values_BER); grid;
xlabel("BER"); ylabel("Sensitivity (dBm)"); title("Sensitivity vs. BER Boundary");
legend(sprintf('Data Rate (kbps) = %i', BW_BB / 1000))

% --- Plot bounding sensitivity vs. data rate ---
Data_rate_values = linspace(1, 50000, res);

num_coeff = 8000 * Vt * sqrt(5 * k * T);
num_ln = log(p_0 * (-2 + 4 * sqrt(0.25 + (BER / (2 * p_0)))));
num = num_coeff * num_ln;

den_coeff = Av * RS * sqrt(CC) * ((1 / n) - 0.5);
den_up_contr = 1 + (8 ./ (2.2 * RD * CC * Data_rate_values));
den_up_contr = (den_up_contr) .^ (3 / 4);
den = den_coeff .* den_up_contr;

P_RF_values_drate = abs(num ./ den);
P_RF_values_drate = 10 * log10(P_RF_values_drate);

subplot(2, 2, 3)
plot(Data_rate_values / 10 ^ 3, P_RF_values_drate); grid;
xlabel("Data Rate (kbps)"); ylabel("Sensitivity (dBm)"); title("Sensitivity vs. Data Rate Boundary");
legend(sprintf('BER = %i', BER))

% --- Plot all three as a surface ---
[DR, BER] = meshgrid(Data_rate_values, BER_values);

num_coeff = 8000 * Vt * sqrt(5 * k * T);
num_ln = log(p_0 * (-2 + 4 * sqrt(0.25 + (BER ./ (2 * p_0)))));
num = num_coeff .* num_ln;

den_coeff = Av * RS * sqrt(CC) * ((1 / n) - 0.5);
den_up_contr = 1 + (8 ./ (2.2 * RD * CC .* DR));
den_up_contr = den_up_contr .^ (3 / 4);
den = den_coeff .* den_up_contr;

P_RF = abs(num ./ den);
P_RF = 10 * log10(P_RF);

subplot(2, 2, 4)
surf(DR / 1e3, BER, P_RF)
shading interp
colorbar

set(gca, 'YScale', 'log')

xlabel("Data Rate (kbps)"); ylabel("BER"); zlabel("Sensitivity (dBm)")
title("Optimal Data Rate / BER / Sensitivity Surface")

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
RS = 50; % Source / antenna resistance (Ohms)
BER = 1e-3; % Desired BER
BW_BB = 1e3; % Baseband signal BW in Hz = Data rate

% Passive voltage gain from matching
Av = 28; % In dB

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
P_RF = 10 ^ (P_RF / 10) * (10 ^ -3);

% Compute the upper bound for the number of allowed stages in the ED
r = 8 / (2.2 * RD * CC * BW_BB);
discriminant = sqrt(1 + r);
N_UPPER_BOUND = 0.25 * discriminant;

% Compute the lower bound for the number of allowed stages in the ED
num_lower_coeff = 5 * k * T;
Q_1_lower = -2 + 4 * ((0.25 + 0.5 * (BER * (((pi + 2) / pi) ^ 0.5) * exp(1 / (pi + 2)))) ^ 0.5);
num_ln = log(Q_1_lower * ((pi / (pi + 2)) ^ 0.5) * exp(-1 / (pi + 2)));
num_ln = num_ln * num_ln;
num_lower = num_lower_coeff * num_ln * Vt;
den_lower = (((1 / n) - 0.5) ^ 2) * (Av * P_RF * RS) ^ 2 * CC;
N_LOWER_BOUND = num_lower / den_lower;
N_LOWER_BOUND = N_LOWER_BOUND ^ (1 / 3);

% Compute upper and lower bounds for ED voltage as a function of N
n_limits = 1:1:2; % Subthreshold slope constant is between 1 and 2
k_limits = 1 ./ (2 .* n_limits * Vt); % Min and Max conversion gain
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

% Select the optimum N
N_opt = (N_LOWER_BOUND + N_UPPER_BOUND) / 2;
fprintf("The average number of stages to meet the given detector criteria is: %0i\n\n", floor(N_opt))

% Summarize design criteria
fprintf("Design criteria summary:\n")
fprintf("Expected RF Input Power (W): %0i\n", P_RF)
fprintf("Desired Bit Rate (Hz): %0i\n", BW_BB)
fprintf("Desired Bit error rate (%%): %0f\n", BER * 100)

% -------------------------------------------------------------------------

% Plot allowed data rate vs. BER for the detector, since both are bounded
% with respect to each other. If this bound is broken in the ED design
% process, the leftmost and rightmist points in the capacity region will
% slide past each other, indicating that a capacity region cannot exist.

% Compute data rate vs. BER bound
BER_values = logspace(log10(0.00001), log10(0.1), 100);
Q_1_lower_sweep = -2 + 4 * ((0.25 + 0.5 * (BER_values * (((pi + 2) / pi) ^ 0.5) * exp(1 / (pi + 2)))) .^ 0.5);
num_ln_sweep = log(Q_1_lower_sweep * ((pi / (pi + 2)) ^ 0.5) * exp(-1 / (pi + 2)));
num_ln_sweep = num_ln_sweep .* num_ln_sweep;
num_lower_sweep = num_lower_coeff .* num_ln_sweep * Vt;
den_lower_sweep = (((1 / n) - 0.5) ^ 2) * (Av * P_RF * RS) ^ 2 * CC;
phi = num_lower_sweep / den_lower_sweep;
BER_function = 16 * ((phi) .^ (2 / 3)) - 1;
DataRate_values = 8 ./ (2.2 * RD * CC * BER_function);

% Plot Data Rate vs BER on log-log scale
figure;
hold on;

% Shade region where data rate <= bound
X = [BER_values * 100, fliplr(BER_values * 100)];
Y = [DataRate_values, ones(size(DataRate_values))*100]; % lower bound (100 bps)
fill(X, Y, [0.8 0.9 1], 'EdgeColor', 'none'); % light blue region

loglog(BER_values * 100, DataRate_values, 'b-', 'LineWidth', 2);
grid on;
xlabel('BER (%)');
ylabel('Data Rate (bps)');
title('Data Rate vs. BER');
xlim([0.001 10]);
ylim([100 1e6]);
set(gca, 'XScale', 'log', 'YScale', 'log');

legend('Feasible Region (≤ bound)', 'Data Rate Upper Bound', 'Location', 'northwest');
hold off;

% --- 3D mapping of sensitivity vs. data rate vs. BER ---

% Define RF power sweep range (in dBm)
P_RF_dBm_values = -80:1:-40; % Sweep from -80 to -40 dBm

% Define BER range (limited)
BER_values = logspace(-4, 0, 100); % 10^-4 to 10^0

% Preallocate for speed
DataRate_matrix = zeros(length(P_RF_dBm_values), length(BER_values));

% Compute data rate for each P_RF value
for i = 1:length(P_RF_dBm_values)
    % Convert dBm to Watts
    P_RF_W = 10^(P_RF_dBm_values(i)/10) * 1e-3;
    Q_1_lower_sweep = -2 + 4 * ((0.25 + 0.5 * (BER_values * (((pi + 2) / pi) ^ 0.5) * exp(1 / (pi + 2)))) .^ 0.5);
    num_ln_sweep = log(Q_1_lower_sweep * ((pi / (pi + 2)) ^ 0.5) * exp(-1 / (pi + 2)));
    num_ln_sweep = num_ln_sweep .* num_ln_sweep;
    num_lower_sweep = num_lower_coeff .* num_ln_sweep * Vt;
    den_lower_sweep = (((1 / n) - 0.5) ^ 2) * (Av * P_RF_W * RS) ^ 2 * CC;
    phi = num_lower_sweep ./ den_lower_sweep;
    BER_function = 16 * ((phi) .^ (2 / 3)) - 1;
    DataRate_values = 8 ./ (2.2 * RD * CC * BER_function);

    % Limit DataRate to 10^1 – 10^6 range
    DataRate_values(DataRate_values < 1e1) = 1e1;
    DataRate_values(DataRate_values > 1e6) = 1e6;

    % Store results
    DataRate_matrix(i, :) = DataRate_values;
end

% Create meshgrid for plotting
[BER_mesh, P_RF_mesh] = meshgrid(BER_values, P_RF_dBm_values);

% 3D plot (BER vs DataRate vs P_RF)
figure;
surf(BER_mesh, DataRate_matrix, P_RF_mesh); % BER (fractional), P_RF (dBm)
xlabel('BER');
ylabel('Data Rate (bps)');
zlabel('P_{RF} (dBm)');
title('3D Relationship: BER vs Data Rate vs RF Input Power');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([1e-4 1e0]);
ylim([1e1 1e6]);
colorbar;
grid on;
view(45,30);
shading interp;

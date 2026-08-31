% ====== Strouhal Number Calculation from Lift Force Data
clear; clc; close all;

% ====== 1. File Selection
% % --- Re 140 --- [un-comment these 2 lines for Re = 140 case]
% U = 0.2169;       % Free-stream velocity [m/s] for Re 140
% filename = '4. re 140 STAR\report_COUPLED_v6\edited\lift_Monitor.csv';
% --- Re 5000 --- [un-comment these 2 lines for Re = 5000 case]
U = 7.844;          % Free-stream velocity [m/s] for Re 5000
filename = '5. re 5000 STAR\report_v8\edited\lift_Monitor.csv';

% ====== 2. Constants
D = 0.01;           % Cylinder diameter [m]
rho = 1.1766;       % Freestream density
q = 0.5 * rho * (U)^2 * D;  % Normalizing factor

% ====== 3. Time mask
% Define the exact time the flow becomes fully developed. This is needed due to sometimes the resulting dominant frequency captured
% is not the true frequency due to the initial transient phase.
t_eval_start = 0.6; 


% ====== 4. Load reference data
raw = readmatrix(filename);
time = raw(:, 1);
lift = raw(:, 2);
cl = lift / q;

% ====== 5. Isolate Stable Oscillation Region
stable_idx = find(time >= t_eval_start, 1, 'first');
    if isempty(stable_idx)
        error('The specified evaluation start time (t=%.2f) exceeds the dataset.', t_eval_start);
    end
% -- Isolate the steady-state portion
t_steady = time(stable_idx:end);
lift_steady = cl(stable_idx:end);
% -- Subtract the mean to remove the offset
% crucial for accurate FFT amplitude
lift_fluct = lift_steady - mean(lift_steady);

% ====== 6. Perform Fast Fourier Transform (FFT)
% Calculate sampling frequency based on actual time steps
dt = mean(diff(t_steady)); 
Fs = 1 / dt;
L = length(t_steady);

% Execute FFT
Y = fft(lift_fluct);

% Compute two-sided spectrum and then single-sided spectrum
P2 = abs(Y / L);
P1 = P2(1:floor(L/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);

% Define the frequency domain
f_domain = Fs * (0:(floor(L/2))) / L;

% Find the dominant frequency (highest amplitude peak)
[~, max_idx] = max(P1);
f_lift = f_domain(max_idx);

% ====== 7. Calculate Strouhal Number
% Strouhal Number Formula: St = (f * D) / U
St = (f_lift * D) / U;

% based on Nelson CC, Nichols RH., 1198, "Comparison of Hybrid Turbulence Models
% for a Circular Cylinder and a Cavity."
t_URANS = (1/f_lift) / 100;
t_DES = (1/f_lift) / 1000;

% ====== 8. Plotting
fprintf('--- Simulation Results ---\n');
fprintf('Lift Oscillation Frequency:   %.3f Hz\n', f_lift);
fprintf('Calculated Strouhal Number:   %.3f\n', St);
fprintf('Suggested URANS timestep:     %.4e\n', t_URANS);
fprintf('Suggested DES timestep:       %.4e\n', t_DES);

% Plot the Lift force history
figure('Name', 'Lift Coefficient Stability Detection', 'Color', 'w');
plot(time, cl, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Lift Coefficient'); hold on;
plot(t_steady, lift_steady, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Analyzed Steady State');
xline(time(stable_idx), 'r--', 'LineWidth', 2, 'Label', sprintf('Evaluation Start (t=%.1fs)', t_eval_start), 'HandleVisibility', 'off');
xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 14);
ylabel('$C_l$', 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter','latex');
title('Lift Coefficient History', 'FontWeight', 'bold', 'FontSize', 14);
legend('Location', 'best');
grid on; box on;

% Plot the Frequency Spectrum
figure('Name', 'FFT Frequency Spectrum', 'Color', 'w');
plot(f_domain, P1, 'b-', 'LineWidth', 1.5);
xlim([0, f_lift*3]); % Zoom in on the relevant frequency range to filter out high-frequency noise
xlabel('Frequency (Hz)', 'FontWeight', 'bold');
ylabel('Amplitude', 'FontWeight', 'bold');
title('FFT Spectrum of Lift Coefficient', 'FontWeight', 'bold');
xline(f_lift, 'r--', sprintf('f_{lift} = %.2f Hz', f_lift), 'LabelVerticalAlignment', 'middle');
grid on; box on;
clear; clc; close all;

% folder_name = '4. re 140 STAR\report_COUPLED_v6\edited';  % BALANCED
folder_name = '5. re 5000 STAR\report_v9\edited';  % BALANCED
num_cvs = 7;    % Process all 7 Control Volumes

% === based on Re ===
% % --- re 140 ---
% Vinf = 0.2196;        % Freestream velocity (m/s)
% f_shed = 8;           % True Vortex Shedding Frequency (Hz)
% t_eval_start = 6;     % Physical time (s) where stable oscillation begins
% --- re 5000 ---
Vinf = 7.844;  % Freestream velocity (m/s)
f_shed = 196;         % True Vortex Shedding Frequency (Hz)
t_eval_start = 0.6;   % Physical time (s) where stable oscillation begins


% 1. Configuration, Normalization, and Styling
n_step = 1;     % Only process every n-th timestep (Set to 1 to read all data)
% --- for energy normalization ---
rho = 1.1766;   % Freestream density
area = 0.01;    % reference area
q = 0.5 * rho * (Vinf)^2 * area;
% --- for time normalization ---
T_shed = 1 / f_shed;  % Shedding Period (s)
num_cycles = 5;       % Exactly 5 full signals to display

% Exact Colors from Final_Exergy_Alex.py
C.Baseline = 'k';                     
C.Sum      = [1 0.647 0];             
C.SumK     = [0.5 0 0.5];             
C.Visc     = 'r';                     
C.Mech     = [0 0.5 0];               
C.Therm    = 'c';                     
C.Vol      = 'b';                     
C.Anablat  = 'm';                     

W.line = 1.6;   % Standardized line width for all plots

% 2. Helper Functions
function [time, val] = get_data(name, col, ref_time, folder, n_skip)
    filename = fullfile(folder, [name, '_Monitor.csv']);
    if exist(filename, 'file')
        data = readmatrix(filename);
        [~,data_column] = size(data);
        if data_column > 2 
            time = data(:,1);
            val = data(:,col+1); 
        else
            time = data(:,1);
            val = data(:,2); 
        end
        time = time(1:n_skip:end);
        val = val(1:n_skip:end);
    else
        warning(['File missing: ', filename, '. Defaulting to 0.']);
        time = ref_time;
        val = zeros(size(ref_time));
    end
end

% 3. Unified Data Loading & Pre-Calculations (All CVs)
[time, drag_baseline] = get_data('drag', 1, [0; 1], folder_name, n_step); 
time_norm = (time - t_eval_start) / T_shed;
cd_baseline = drag_baseline / q;
eval_mask = (time_norm >= 0) & (time_norm <= num_cycles);
mean_baseline = mean(cd_baseline(eval_mask), 'omitnan');

% Preallocate cell arrays for calculated coefficients
Cd_phi = cell(1, num_cvs); Cd_Etp = cell(1, num_cvs); Cd_theta = cell(1, num_cvs); 
Cd_dek = cell(1, num_cvs); Cd_Energy = cell(1, num_cvs);
Cd_aphi = cell(1, num_cvs); Cd_xm = cell(1, num_cvs); Cd_xth = cell(1, num_cvs);
Cd_anabla = cell(1, num_cvs); Cd_Ex_Store = cell(1, num_cvs); Cd_Ex_Kinetic = cell(1, num_cvs);
Cd_Exergy = cell(1, num_cvs);
Cd_dEnth = cell(1, num_cvs); Cd_EnthFlux_h = cell(1, num_cvs); 
Cd_EnthFlux_p = cell(1, num_cvs); Cd_EnthFlux_tau = cell(1, num_cvs);
Cd_Enthalpy = cell(1, num_cvs);


% Preallocate for Energy Delta Plot
cv_temp_aphi = cell(1, num_cvs);
cv_temp_xvt = cell(1, num_cvs);
cv_temp_Xm = cell(1, num_cvs);
avg_XKE_CV = zeros(1, num_cvs);
avg_Aphi_CV = zeros(1, num_cvs);

% Preallocate Error Arrays
err_Energy = zeros(1, num_cvs);
err_Exergy = zeros(1, num_cvs);
err_Enthalpy = zeros(1, num_cvs);

for c = 1:num_cvs
    % Direct non-cumulative data loading for each specific CV column/file
    load_data = @(name) get_data(name, c, time, folder_name, n_step);
    
    % Volumetric Terms (Non-cumulative)
    [~, phi_vol]    = load_data('phi_vol');
    [~, theta_vol]  = load_data('theta_vol');
    [~, ek_vol]     = load_data('ek_vol');
    [~, aphi_vol]   = load_data('aphi_vol');
    [~, anabla_vol] = load_data('anabla_vol');
    [~, massspec]   = load_data('massspec_vol');
    [~, enth_storage_vol] = load_data('enthalpy_storage_vol'); % NEW: Total enthalpy volume storage
    [~, phi_alt_vol] = load_data('phi_alt_vol');   % NEW: Phi volume term based on internal energy
    [~, aphi_alt_vol] = load_data('aphi_alt_vol');   % NEW: Phi volume term based on internal energy

    dek_dt = gradient(ek_vol, time);
    dxv_dt = gradient(massspec, time);
    d_enth_dt = gradient(enth_storage_vol, time);
    d_phi_alt = gradient(phi_alt_vol, time);
    d_aphi_alt = gradient(aphi_alt_vol, time);

    % Mechanical Fluxes (Non-cumulative)
    [~, EA_bot] = load_data('ea_bot'); [~, EA_top] = load_data('ea_top');
    [~, EA_left] = load_data('ea_left'); [~, EA_right] = load_data('ea_right');
    EA_total = EA_bot + EA_top + EA_left + EA_right;

    [~, EV_bot] = load_data('ev_bot'); [~, EV_top] = load_data('ev_top');
    [~, EV_left] = load_data('ev_left'); [~, EV_right] = load_data('ev_right');
    EV_total = EV_bot + EV_top + EV_left + EV_right;

    [~, EP_bot] = load_data('ep_bot'); [~, EP_top] = load_data('ep_top');
    [~, EP_left] = load_data('ep_left'); [~, EP_right] = load_data('ep_right');
    EP_total = EP_bot + EP_top + EP_left + EP_right;

    [~, Etau_bot] = load_data('etau_bot'); [~, Etau_top] = load_data('etau_top');
    [~, Etau_left] = load_data('etau_left'); [~, Etau_right] = load_data('etau_right');
    Etau_total = Etau_bot + Etau_top + Etau_left + Etau_right;

    [~, Etp_bot] = load_data('etp_bot'); [~, Etp_top] = load_data('etp_top');
    [~, Etp_left] = load_data('etp_left'); [~, Etp_right] = load_data('etp_right');
    Etp_total = Etp_bot + Etp_left + Etp_top + Etp_right;
    % Etp_total = EA_total + EV_total + EP_total - Etau_total;
    
    % Total Enthalpy Wake Flux Components (Method 3)
    [~, Hflux_bot]   = load_data('enth_flux_enthalpy_bot'); [~, Hflux_top]   = load_data('enth_flux_enthalpy_top');
    [~, Hflux_left]  = load_data('enth_flux_enthalpy_left'); [~, Hflux_right] = load_data('enth_flux_enthalpy_right');
    Hflux_total      = Hflux_bot + Hflux_top + Hflux_left + Hflux_right;

    [~, Pflux_bot]   = load_data('enth_flux_pressure_bot'); [~, Pflux_top]   = load_data('enth_flux_pressure_top');
    [~, Pflux_left]  = load_data('enth_flux_pressure_left'); [~, Pflux_right] = load_data('enth_flux_pressure_right');
    Pflux_total      = Pflux_bot + Pflux_top + Pflux_left + Pflux_right;

    [~, Tauflux_bot]   = load_data('enth_flux_viscous_bot'); [~, Tauflux_top]   = load_data('enth_flux_viscous_top');
    [~, Tauflux_left]  = load_data('enth_flux_viscous_left'); [~, Tauflux_right] = load_data('enth_flux_viscous_right');
    Tauflux_total      = Tauflux_bot + Tauflux_top + Tauflux_left + Tauflux_right;

    % Thermal Exergy Flux (Non-cumulative)
    [~, Xth_bot] = load_data('xth_bot'); [~, Xth_top] = load_data('xth_top');
    [~, Xth_left] = load_data('xth_left'); [~, Xth_right] = load_data('xth_right');
    Xth_total = Xth_bot + Xth_top + Xth_left + Xth_right;
    
    % --- New phi flux !!! ---
    [~, phi_alt_bot] = load_data('phi_alt_flux_bot'); [~, phi_alt_top] = load_data('phi_alt_flux_top');
    [~, phi_alt_left] = load_data('phi_alt_flux_left'); [~, phi_alt_right] = load_data('phi_alt_flux_right');
    phi_alt_total = phi_alt_bot + phi_alt_top + phi_alt_left + phi_alt_right;
    phi_alt = d_phi_alt + phi_alt_total;
    % --- New Aphi flux !!! ---
    [~, aphi_alt_bot] = load_data('aphi_alt_flux_bot'); [~, aphi_alt_top] = load_data('aphi_alt_flux_top');
    [~, aphi_alt_left] = load_data('aphi_alt_flux_left'); [~, aphi_alt_right] = load_data('aphi_alt_flux_right');
    aphi_alt_total = aphi_alt_bot + aphi_alt_top + aphi_alt_left + aphi_alt_right;
    aphi_alt = d_aphi_alt + aphi_alt_total;

    
    % Data Logging for Delta Transfers
    cv_temp_aphi{c} = aphi_vol;
    % cv_temp_aphi{c} = aphi_alt;         % with new Aphi formulation
    cv_temp_xvt{c}  = dek_dt + dxv_dt;
    cv_temp_Xm{c}   = Etp_total;

    % --- Method 1: Energy ---
    Cd_phi{c} = (phi_vol / Vinf) / q;     % default formulation
    % Cd_phi{c} = (phi_alt / Vinf) / q;       % with New Phi formulation
    Cd_Etp{c} = (Etp_total / Vinf) / q;
    Cd_theta{c} = (theta_vol / Vinf) / q;
    Cd_dek{c} = (dek_dt / Vinf) / q;
    Cd_Energy{c} = ((dek_dt + Etp_total + theta_vol + phi_vol) / Vinf) / q;
    
    mean_energy = mean(Cd_Energy{c}(eval_mask), 'omitnan');
    err_Energy(c) = abs(mean_energy - mean_baseline) / mean_baseline * 100;

    % --- Method 2: Exergy ---
    Exergy_Storage = dxv_dt + dek_dt;
    Exergy_Mech_Flux = Etp_total;
    Exergy_Boundary_Flux = Xth_total + Exergy_Mech_Flux;
    Exergy_Kinetic = Exergy_Mech_Flux + Exergy_Storage;

    Cd_aphi{c} = (aphi_vol / Vinf) / q;
    Cd_anabla{c} = (anabla_vol / Vinf) / q;
    Cd_xth{c} = (Xth_total / Vinf) / q;
    Cd_xm{c} = (Exergy_Mech_Flux / Vinf) / q;
    Cd_Ex_Store{c} = (Exergy_Storage / Vinf) / q;
    Cd_Ex_Kinetic{c} = (Exergy_Kinetic / Vinf) / q;
    Cd_Exergy{c} = ((Exergy_Storage + Exergy_Boundary_Flux + aphi_vol + anabla_vol) / Vinf) / q;
    
    mean_exergy = mean(Cd_Exergy{c}(eval_mask), 'omitnan');
    err_Exergy(c) = abs(mean_exergy - mean_baseline) / mean_baseline * 100;
    

    % --- Method 3: Enthalpy Balance (Total Energy) ---
    Cd_dEnth{c}        = (d_enth_dt / Vinf) / q;
    Cd_EnthFlux_h{c}   = (Hflux_total / Vinf) / q;
    Cd_EnthFlux_p{c}   = (Pflux_total / Vinf) / q;
    Cd_EnthFlux_tau{c} = (Tauflux_total / Vinf) / q;
    Cd_Enthalpy{c}     = ((d_enth_dt + Hflux_total + Pflux_total + Tauflux_total) / Vinf) / q;

    mean_enthalpy = mean(Cd_Enthalpy{c}(eval_mask), 'omitnan');
    err_Enthalpy(c) = abs(mean_enthalpy - mean_baseline) / mean_baseline * 100;

    % Delta Means Processing
    temp_Cd_XKE  = ((cv_temp_Xm{c} + cv_temp_xvt{c}) / Vinf) / q;
    temp_Cd_Aphi = (cv_temp_aphi{c} / Vinf) / q;
    avg_XKE_CV(c)  = mean(temp_Cd_XKE(eval_mask), 'omitnan');
    avg_Aphi_CV(c) = mean(temp_Cd_Aphi(eval_mask), 'omitnan');
end

% PLOTTING SECTIONS (TABULATED)
%% 1. Plot Method 1: Energy Balance
fig1 = figure('Name', 'Method 1: Energy Balance', 'Position', [100, 100, 750, 520], 'Color', 'w');
tg1 = uitabgroup(fig1);

for c = 1:num_cvs
    t = uitab(tg1, 'Title', sprintf('CV %d', c));
    ax = axes('Parent', t); hold(ax, 'on');

    plot(ax, time_norm, Cd_phi{c}, 'Color', C.Visc, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_\Phi$');
    plot(ax, time_norm, Cd_Etp{c}, 'Color', C.Mech, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{E}_{TP}}$');
    plot(ax, time_norm, Cd_theta{c}, 'Color', C.Therm, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_\Theta$');
    plot(ax, time_norm, Cd_dek{c}, 'Color', C.Vol, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{dE_k/dt}$');
    plot(ax, time_norm, Cd_Energy{c}, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_D$');
    plot(ax, time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');

    title(ax, sprintf('Energy Balance Drag Decomposition -- CV%d (Error: %.2f%%)', c, err_Energy(c)));
    xlabel(ax, '$t/T$', 'Interpreter', 'latex', 'FontSize', 12); 
    ylabel(ax, 'Energy Coefficient', 'FontSize', 12);
    xlim(ax, [0, num_cycles]); 
    grid(ax, 'on'); ax.GridAlpha = 0.3; ax.XTick = 0:num_cycles;
    legend(ax, 'Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 14, 'Box', 'off');
end

%% 2. Plot Method 2: Exergy Balance
fig2 = figure('Name', 'Method 2: Exergy Balance', 'Position', [200, 100, 750, 520], 'Color', 'w');
tg2 = uitabgroup(fig2);

for c = 1:num_cvs
    t = uitab(tg2, 'Title', sprintf('CV %d', c));
    ax = axes('Parent', t); hold(ax, 'on');

    plot(ax, time_norm, Cd_aphi{c}, 'Color', C.Visc, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{A}_\Phi}$');
    plot(ax, time_norm, Cd_xm{c}, 'Color', C.Mech, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{X}_m}$');
    plot(ax, time_norm, Cd_xth{c}, 'Color', C.Therm, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{X}_{th}}$');
    plot(ax, time_norm, Cd_anabla{c}, 'Color', C.Anablat, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{A_{\nabla T}}$');
    plot(ax, time_norm, Cd_Ex_Store{c}, 'Color', C.Vol, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\dot{X}_{v(t)}}$');
    plot(ax, time_norm, Cd_Ex_Kinetic{c}, 'Color', C.SumK, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', '$\Sigma X_k$');
    plot(ax, time_norm, Cd_Exergy{c}, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_{\dot{X}}$');
    plot(ax, time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');

    title(ax, sprintf('Far-Field Exergy Balance -- CV%d (Error: %.2f%%)', c, err_Exergy(c)));
    xlabel(ax, '$t/T$', 'Interpreter', 'latex', 'FontSize', 16); 
    ylabel(ax, 'Exergy Coefficient', 'FontSize', 16);
    xlim(ax, [0, num_cycles]);
    grid(ax, 'on'); ax.GridAlpha = 0.3; ax.XTick = 0:num_cycles;
    legend(ax, 'Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 16, 'Box', 'off');
end

%% 3. Plot Method 3: Enthalpy Balance
fig3 = figure('Name', 'Method 3: Enthalpy Balance', 'Position', [300, 100, 750, 520], 'Color', 'w');
tg3 = uitabgroup(fig3);

for c = 1:num_cvs
    t = uitab(tg3, 'Title', sprintf('CV %d', c));
    ax = axes('Parent', t); hold(ax, 'on');

    plot(ax, time_norm, Cd_dEnth{c}, 'Color', C.Vol, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{dH_0/dt}$');
    plot(ax, time_norm, Cd_EnthFlux_h{c}, 'Color', C.Therm, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\Delta\dot{H}_0}$');
    plot(ax, time_norm, Cd_EnthFlux_p{c}, 'Color', C.Mech, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\Delta\dot{W}_p}$');
    plot(ax, time_norm, Cd_EnthFlux_tau{c}, 'Color', C.Visc, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$C_{\Delta\dot{W}_\tau}$');
    plot(ax, time_norm, Cd_Enthalpy{c}, 'Color', C.Sum, 'LineStyle', '-', 'LineWidth', W.line, 'DisplayName', '$\Sigma C_D$');
    plot(ax, time_norm, cd_baseline, 'Color', C.Baseline, 'LineStyle', '--', 'LineWidth', W.line, 'DisplayName', 'Near Field $C_d$');

    title(ax, sprintf('Enthalpy Balance Drag Decomposition -- CV%d (Error: %.2f%%)', c, err_Enthalpy(c)));
    xlabel(ax, '$t/T$', 'Interpreter', 'latex', 'FontSize', 16); 
    ylabel(ax, 'Enthalpy Coefficient', 'FontSize', 16);
    xlim(ax, [0, num_cycles]); 
    grid(ax, 'on'); ax.GridAlpha = 0.3; ax.XTick = 0:num_cycles;
    legend(ax, 'Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 16, 'Box', 'off');
end

%% 4. ENERGY DELTAS (BETWEEN CONTROL VOLUMES)
dKE   = zeros(1, num_cvs - 1);
dAphi = zeros(1, num_cvs - 1);
label_names = cell(1, num_cvs - 1);

for i = 1:(num_cvs - 1)
    dKE(i) = avg_XKE_CV(i) - avg_XKE_CV(i+1);
    dAphi(i) = avg_Aphi_CV(i+1) - avg_Aphi_CV(i);
    label_names{i} = sprintf('CV%d \\rightarrow CV%d', i, i+1);
end

labels = categorical(label_names);
labels = reordercats(labels, label_names);

figure("Name", "Energy Transfer Imbalance", 'Position', [500, 100, 900, 520], 'Color', 'w');
hold on; box on;

b = bar(labels, [dKE' dAphi'], 'LineWidth', 1.0, 'BarWidth', 1);

b(1).FaceColor = [0.2 0.2 0.2];     
b(2).FaceColor = [0.8 0.2 0.2];     
b(1).EdgeColor = 'k';
b(2).EdgeColor = 'k';

ylabel('$\Delta C_{[]}$', 'Interpreter', 'latex', 'FontSize', 14);
title('Energy transfer between CV', 'FontWeight', 'bold');
ax = gca;
ax.FontSize = 14; 
grid on; ax.GridAlpha = 0.3;

legend({'$\Delta X_{KE}$ lost', '$\Delta A_{\Phi}$ gained'}, ...
       'Interpreter', 'latex', 'FontSize', 11, 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

xtips1 = b(1).XEndPoints; xtips2 = b(2).XEndPoints;
ytips1 = b(1).YEndPoints; ytips2 = b(2).YEndPoints;

for i = 1:length(dKE)
    pct = (dAphi(i) - dKE(i)) / abs(dKE(i)) * 100;
    x_center = (xtips1(i) + xtips2(i)) / 2;
    y_max = max(ytips1(i), ytips2(i));
    
    if isnan(pct)
        label_text = 'n/a';
    else
        label_text = sprintf('%+.1f%%', pct);
    end
    
    text(x_center, y_max, label_text, 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
end

ylim([0, max(max(ytips1, ytips2)) * 1.15]);
% Exergy Comparison: STAR-CCM+ vs. ANSYS Fluent
clear; clc; close all;

% % === Re 140 case ===
% folder_star   = '4. re 140 STAR\report_COUPLED_v5.5C_thalf_constRho\edited';
% folder_fluent = '4. re 140 STAR\report_FLUENT_COUPLED_ALEX\data\edited';
% f_shed = 8;       % --- shedding frequency ---
% Vinf = 0.2196;    % --- freestream velocity ---
% % --- averaging time window ---
% t_start = 6.0;
% t_end   = 7.0;


% === Re 5000 case ===
folder_star   = '5. re 5000 STAR\report_v8\edited';
folder_fluent = '6. re 5000 FLUENT\data\edited';
f_shed   = 197;    % --- shedding frequency ---
Vinf     = 7.844;  % --- freestream velocity ---
% --- averaging time window ---
t_start = 0.6;
t_end   = 1.0;


n_cv = 7;
% -- Cycle-normalization for the Exergy Balance tab plots
T_shed       = 1 / f_shed;     % Shedding Period (s)
t_eval_start = t_start;        % Physical time (s) where stable oscillation begins
num_cycles   = 5;              % Exactly 5 full signals to display


rho = 1.1766;
area = 0.01;
q = 0.5 * rho * (Vinf)^2 * area;
power_ref = q * Vinf;

terms = {'Aphi', 'Anabla', 'Xm', 'Xth', 'Xv', 'Xk', 'Total'};
term_labels = {'Viscous Anergy', 'Thermal Anergy', 'Mech Exergy', 'Thermal Exergy', 'Unsteady Storage', 'Kinetic Exergy', 'Total'};

% Exact Colors from Final_Exergy_Alex.py (matching drag_compare_STAR_4_fix_step_ALEX.m)
C.Baseline = 'k';                     % Black
C.Sum      = [1 0.647 0];             % Orange
C.SumK     = [0.5 0 0.5];             % Purple
C.Visc     = 'r';                     % Red (aphi)
C.Mech     = [0 0.5 0];               % Dark Green (xm) matching Fluent #008000
C.Therm    = 'c';                     % Cyan (xth)
C.Vol      = 'b';                     % Blue (dXv)
C.Anablat  = 'm';                     % Magenta (anablat)

W.line = 1.6;   % Standardized line width for all plots

% Order matches fields = {'aphi','anab','xm','xth','xv','xk','tot'} used in plot_decomp_tab
term_colors = {C.Visc, C.Anablat, C.Mech, C.Therm, C.Vol, C.SumK, C.Sum};
term_linestyles = {'-', '-', '-', '-', '-', '--', '-'};   % Xk dashed, matching Cd_Ex_Kinetic in the reference script
term_labels_latex = {'$\dot{\mathcal{A}}_\Phi$', '$\mathcal{A}_{\nabla T}$', '$\dot{\mathcal{X}}_m$', ...
                     '$\dot{\mathcal{X}}_{th}$', '$\dot{\mathcal{X}}_{v(t)}$', '$\Sigma \mathcal{X}_k$', '$\Sigma \dot{\mathcal{X}}$'};
% Same symbols, with the leading "C" (coefficient) dropped -- used only for
% the bar chart legend below, not the decomposition line plots.
bar_labels_latex = {'$\dot{\mathcal{A}}_\Phi$', '$\mathcal{A}_{\nabla T}$', '$\dot{\mathcal{X}}_m$', ...
                     '$\dot{\mathcal{X}}_{th}$', '$\dot{\mathcal{X}}_{v(t)}$', '$\Sigma \mathcal{X}_k$', '$\Sigma \dot{\mathcal{X}}$'};

star_matrix   = zeros(length(terms), n_cv);
fluent_matrix = zeros(length(terms), n_cv);
S_plot = cell(n_cv, 1); F_plot = cell(n_cv, 1);

% Arrays to store the exact phase-aligned start times for each solver per CV
t_start_S = zeros(n_cv, 1); 
t_start_F = zeros(n_cv, 1);

for cv = 1:n_cv
    fprintf('--- Processing CV %d ---\n', cv);

    % --- Star-CCM+: load + compute using the STAR-specific formulation ---
    S = compute_star_exergy_cv(cv, folder_star, power_ref);

    % --- Fluent: load + compute using the FLUENT-specific formulation ---
    F = compute_fluent_exergy_cv(cv, folder_fluent, power_ref);

    % --- Phase Alignment for STAR-CCM+ ---
    idx_s = find(S.t >= t_start);
    tot_s_ac = S.tot(idx_s) - mean(S.tot(idx_s), 'omitnan');
    % Find the first positive-going zero crossing
    cross_s = find(tot_s_ac(1:end-1) <= 0 & tot_s_ac(2:end) > 0, 1);
    t_start_S(cv) = interp1(tot_s_ac(cross_s:cross_s+1), S.t(idx_s(cross_s:cross_s+1)), 0);
    
    % New dynamically shifted mask bounding exactly num_cycles
    mask_s = (S.t >= t_start_S(cv)) & (S.t <= (t_start_S(cv) + num_cycles * T_shed));

    star_matrix(1, cv) = mean(S.aphi(mask_s), 'omitnan');
    star_matrix(2, cv) = mean(S.anab(mask_s), 'omitnan');
    star_matrix(3, cv) = mean(S.xm(mask_s),   'omitnan');
    star_matrix(4, cv) = mean(S.xth(mask_s),  'omitnan');
    star_matrix(5, cv) = mean(S.xv(mask_s),   'omitnan');
    star_matrix(6, cv) = mean(S.xk(mask_s),   'omitnan');
    star_matrix(7, cv) = mean(S.tot(mask_s),  'omitnan');

    % --- Phase Alignment for FLUENT ---
    idx_f = find(F.t >= t_start);
    tot_f_ac = F.tot(idx_f) - mean(F.tot(idx_f), 'omitnan');
    % Find the first positive-going zero crossing
    cross_f = find(tot_f_ac(1:end-1) <= 0 & tot_f_ac(2:end) > 0, 1);
    t_start_F(cv) = interp1(tot_f_ac(cross_f:cross_f+1), F.t(idx_f(cross_f:cross_f+1)), 0);
    
    % New dynamically shifted mask bounding exactly num_cycles
    mask_f = (F.t >= t_start_F(cv)) & (F.t <= (t_start_F(cv) + num_cycles * T_shed));

    fluent_matrix(1, cv) = mean(F.aphi(mask_f), 'omitnan');
    fluent_matrix(2, cv) = mean(F.anab(mask_f), 'omitnan');
    fluent_matrix(3, cv) = mean(F.xm(mask_f),   'omitnan');
    fluent_matrix(4, cv) = mean(F.xth(mask_f),  'omitnan');
    fluent_matrix(5, cv) = mean(F.xv(mask_f),   'omitnan');
    fluent_matrix(6, cv) = mean(F.xk(mask_f),   'omitnan');
    fluent_matrix(7, cv) = mean(F.tot(mask_f),  'omitnan');

    % --- Store full time series (already coefficient-normalized) for the
    %     tabbed decomposition plots below ---
    S_plot{cv} = S;
    F_plot{cv} = F;
end

%% ========================================================================
%  TASK 1: RELATIVE DIFFERENCE BAR GRAPH (Star vs Fluent), ALL 7 CVs
%  ========================================================================
% rel_diff(term, cv) = (Star_term - Fluent_term) / |Fluent_Total(cv)| * 100
% Every term (for a given CV) is normalized by that CV's own Fluent Total
% exergy magnitude, rather than by each term's own (sometimes near-zero)
% mean. This removes the near-zero-denominator blow-ups seen with Xth/Xv
% (and now Xm/Xk), but changes what the metric means: it now reads as
% "this term's mismatch, as a % of the CV's overall exergy budget" rather
% than "this term's mismatch, as a % of its own typical size". A term
% that's small relative to Total will always show a small % here even if
% Star and Fluent disagree substantially about that term specifically.
idx_total = find(strcmp(terms, 'Total'));
denom = abs(fluent_matrix(idx_total, :));   % 1 x n_cv, per-CV Fluent Total magnitude

rel_diff = (star_matrix - fluent_matrix) ./ denom * 100;   % denom broadcasts over all 7 term rows

fprintf('Percentage Differences (Star-CCM+ vs Fluent) relative to Fluent Total Exergy\n');
fprintf('====================================================================================================\n');

% Print header: Columns for each CV[cite: 1]
fprintf('%-15s', 'Term \ CV');
for cv = 1:n_cv
    fprintf('| %-10s', sprintf('CV %d', cv));
end
fprintf('\n----------------------------------------------------------------------------------------------------\n');

% Print rows for each term[cite: 1]
for i = 1:length(terms)
    % Add a separator before the last row (Total Exergy)[cite: 1]
    if i == length(terms)
        fprintf('----------------------------------------------------------------------------------------------------\n');
    end

    fprintf('%-15s', terms{i});
    for cv = 1:n_cv
        fprintf('| %-9.4f%%', rel_diff(i, cv));
    end
    fprintf('\n');
end
fprintf('====================================================================================================\n\n');


cv_labels = {'2D', '5D', '10D', '15D', '20D', '30D', '50D'};
cv_cat = categorical(cv_labels);
cv_cat = reordercats(cv_cat, cv_labels);

figure('Name', 'Relative Exergy Difference by CV', 'Position', [100, 100, 1200, 600], 'Color', 'w');
% rel_diff' is (CV x term): one x-axis group per CV, one bar per term within each group
b = bar(cv_cat, rel_diff', 'grouped');

% Match bar colors to the Exergy Balance line-plot palette (term_colors is
% already in the same term order: aphi, anab, xm, xth, xv, xk, tot)
for k = 1:numel(term_colors)
    b(k).FaceColor = term_colors{k};
    b(k).EdgeColor = 'k';
end

ylabel('Difference (% of CV Total Exergy)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Streamwise CV Length', 'FontSize', 14, 'FontWeight', 'bold');
yticks('auto');
title('Exergy Term Mismatch, Normalized by Total Exergy (Star-CCM+ vs Fluent Reference)', 'FontSize', 16);
% subtitle('Each term normalized by that CV''s mean', 'FontSize', 9, 'FontAngle', 'italic');
grid on; ax = gca; ax.GridAlpha = 0.3;
ax.YTick = floor(ax.YLim(1)):5:ceil(ax.YLim(2));

% Legend combines the romanized name with its latex symbol (no "C" prefix)
combined_legend = cell(1, numel(term_labels));
for k = 1:numel(term_labels)
    combined_legend{k} = sprintf('%s  %s', term_labels{k}, bar_labels_latex{k});
end
legend(combined_legend, 'Interpreter', 'latex', 'Location', 'southoutside', ...
       'Orientation', 'horizontal', 'NumColumns', 4, 'Box', 'off', 'FontSize', 12);

%% ========================================================================
%  TASK 2: TABBED DRAG-DECOMPOSITION FIGURES (Star-CCM+ and Fluent), 7 CVs
%  X-axis is cycle-normalized (t/T), showing exactly num_cycles periods,
%  using the same logic as drag_compare_STAR_4_fix_step_ALEX.m.
%  ========================================================================

fig_star = uifigure('Name', 'Exergy Decomposition: Star-CCM+', ...
                     'Position', [100, 80, 950, 680], 'Color', 'w');
tg_star = uitabgroup(fig_star, 'Units', 'normalized', 'Position', [0 0 1 1]);
for cv = 1:n_cv
    tab = uitab(tg_star, 'Title', sprintf('CV%d', cv));
    ax = uiaxes(tab, 'Units', 'normalized', 'Position', [0.08 0.15 0.88 0.75]);
    plot_decomp_tab(ax, S_plot{cv}, term_labels_latex, term_colors, term_linestyles, W, cv, ...
                     'Star-CCM+', t_start_S(cv), T_shed, num_cycles);
end

fig_fluent = uifigure('Name', 'Exergy Decomposition: Fluent', ...
                       'Position', [200, 80, 950, 680], 'Color', 'w');
tg_fluent = uitabgroup(fig_fluent, 'Units', 'normalized', 'Position', [0 0 1 1]);
for cv = 1:n_cv
    tab = uitab(tg_fluent, 'Title', sprintf('CV%d', cv));
    ax = uiaxes(tab, 'Units', 'normalized', 'Position', [0.08 0.15 0.88 0.75]);
    plot_decomp_tab(ax, F_plot{cv}, term_labels_latex, term_colors, term_linestyles, W, cv, ...
                     'Fluent', t_start_F(cv), T_shed, num_cycles);
end


%  LOCAL FUNCTIONS  (must appear after all script-level commands)
% ========================================================================

% -- Data loader (file-exists guard + cv_idx+1 offset, matching the
%    convention already validated in the original drag-decomposition script) --
function [time, val] = load_data_robust(folder, filename, cv_idx)
    filepath = fullfile(folder, [filename, '.csv']);
    if ~exist(filepath, 'file')
        warning(['File missing: ', filename, '. Using zeros.']);
        time = zeros(100, 1); val = zeros(100, 1); return;
    end
    matrix_data = readmatrix(filepath);

    % Exception: If only 2 columns exist (Flow time, Left flux), return column 2
    if size(matrix_data, 2) == 2
        time = matrix_data(:, 1);
        val = matrix_data(:, 2);
    else
        % Standard: Col 1: Time, Col 2: CV1 (Index cv_idx + 1)
        time = matrix_data(:, 1);
        target_col = cv_idx + 1;

        if target_col > size(matrix_data, 2)
            target_col = size(matrix_data, 2);
        end
        val = matrix_data(:, target_col);
    end
end

% -- Cumulative loader for Fluent "_vol" monitors --------------------------
% Fluent's *_vol monitor files store each CV column as an incremental
% (non-cumulative) value, so the true CVn quantity is the running sum of
% raw columns CV1..CVn (e.g. true CV2 = raw CV1 + raw CV2, true CV3 =
% raw CV1 + raw CV2 + raw CV3, ...). This loader reproduces that
% cumulative-sum logic; used only for Fluent's aphi_vol / anabla_vol /
% ek_vol / massspec_vol monitors.
function [time, val] = load_data_robust_cumulative(folder, filename, cv_idx)
    filepath = fullfile(folder, [filename, '.csv']);
    if ~exist(filepath, 'file')
        warning(['File missing: ', filename, '. Using zeros.']);
        time = zeros(100, 1); val = zeros(100, 1); return;
    end
    matrix_data = readmatrix(filepath);

    if size(matrix_data, 2) == 2
        % Only one CV column present -- nothing to accumulate
        time = matrix_data(:, 1);
        val = matrix_data(:, 2);
    else
        time = matrix_data(:, 1);
        last_col = cv_idx + 1;
        if last_col > size(matrix_data, 2)
            last_col = size(matrix_data, 2);
        end
        % Cumulative sum of raw CV1..CVcv_idx columns
        val = sum(matrix_data(:, 2:last_col), 2);
    end
end

% -- STAR-CCM+ Exergy Balance formulation --------------------------------
% Mechanical exergy flux is built from the four separated flux components
% (EA, EV, EP, Etau) on each face, per the Star-CCM+ Separated Method used
% throughout drag_compare_STAR_4_fix_step_ALEX.m.
function S = compute_star_exergy_cv(cv, folder, power_ref)
    [t, aphi]   = load_data_robust(folder, 'aphi_vol_Monitor', cv);
    [~, anab]   = load_data_robust(folder, 'anabla_vol_Monitor', cv);
    [~, ek]     = load_data_robust(folder, 'ek_vol_Monitor', cv);
    [~, mass]   = load_data_robust(folder, 'massspec_vol_Monitor', cv);

    [~, EA_b]   = load_data_robust(folder, 'EA_bot_Monitor', cv);
    [~, EA_t]   = load_data_robust(folder, 'EA_top_Monitor', cv);
    [~, EA_l]   = load_data_robust(folder, 'EA_left_Monitor', cv);
    [~, EA_r]   = load_data_robust(folder, 'EA_right_Monitor', cv);

    [~, EV_b]   = load_data_robust(folder, 'EV_bot_Monitor', cv);
    [~, EV_t]   = load_data_robust(folder, 'EV_top_Monitor', cv);
    [~, EV_l]   = load_data_robust(folder, 'EV_left_Monitor', cv);
    [~, EV_r]   = load_data_robust(folder, 'EV_right_Monitor', cv);

    [~, EP_b]   = load_data_robust(folder, 'EP_bot_Monitor', cv);
    [~, EP_t]   = load_data_robust(folder, 'EP_top_Monitor', cv);
    [~, EP_l]   = load_data_robust(folder, 'EP_left_Monitor', cv);
    [~, EP_r]   = load_data_robust(folder, 'EP_right_Monitor', cv);

    [~, Etau_b] = load_data_robust(folder, 'Etau_bot_Monitor', cv);
    [~, Etau_t] = load_data_robust(folder, 'Etau_top_Monitor', cv);
    [~, Etau_l] = load_data_robust(folder, 'Etau_left_Monitor', cv);
    [~, Etau_r] = load_data_robust(folder, 'Etau_right_Monitor', cv);
    
    [~, Etp_b] = load_data_robust(folder, 'Etp_bot_Monitor', cv);
    [~, Etp_t] = load_data_robust(folder, 'Etp_top_Monitor', cv);
    [~, Etp_l] = load_data_robust(folder, 'Etp_left_Monitor', cv);
    [~, Etp_r] = load_data_robust(folder, 'Etp_right_Monitor', cv);

    [~, Xth_b]  = load_data_robust(folder, 'Xth_bot_Monitor', cv);
    [~, Xth_t]  = load_data_robust(folder, 'Xth_top_Monitor', cv);
    [~, Xth_l]  = load_data_robust(folder, 'Xth_left_Monitor', cv);
    [~, Xth_r]  = load_data_robust(folder, 'Xth_right_Monitor', cv);

    dek_dt = gradient(ek, t);
    dxv_dt = gradient(mass, t);

    % xm  = (EA_b + EA_t + EA_l + EA_r) + (EV_b + EV_t + EV_l + EV_r) + ...
    %       (EP_b + EP_t + EP_l + EP_r) - (Etau_b + Etau_t + Etau_l + Etau_r);
    xm  = Etp_b + Etp_t + Etp_l + Etp_r;
    xv  = dek_dt + dxv_dt;
    xth = Xth_b + Xth_t + Xth_l + Xth_r;
    xk  = xm + xv;
    tot = aphi + anab + xm + xth + xv;

    S.t    = t;
    S.aphi = aphi / power_ref;
    S.anab = anab / power_ref;
    S.xm   = xm   / power_ref;
    S.xth  = xth  / power_ref;
    S.xv   = xv   / power_ref;
    S.xk   = xk   / power_ref;
    S.tot  = tot  / power_ref;
end

% -- Fluent Exergy Balance formulation ------------------------------------
% Fluent's report set exports the mechanical exergy flux already combined
% per face ('etp'), unlike Star's four separated (EA/EV/EP/Etau) terms, so
% xm is built directly from the etp faces rather than re-derived from
% separated components. The "_vol" monitors (aphi/anablat/ek/massspec) are
% stored incrementally per CV, so they must be loaded cumulatively
% (CVn = raw CV1 + ... + raw CVn) via load_data_robust_cumulative.
function F = compute_fluent_exergy_cv(cv, folder, power_ref)
    [t, aphi]   = load_data_robust_cumulative(folder, 'aphi_vol_Monitor', cv);
    [~, anab]   = load_data_robust_cumulative(folder, 'anabla_vol_Monitor', cv);
    [~, ek]     = load_data_robust_cumulative(folder, 'ek_vol_Monitor', cv);
    [~, mass]   = load_data_robust_cumulative(folder, 'massspec_vol_Monitor', cv);

    [~, Etp_b]  = load_data_robust(folder, 'etp_bot_Monitor', cv);
    [~, Etp_t]  = load_data_robust(folder, 'etp_top_Monitor', cv);
    [~, Etp_l]  = load_data_robust(folder, 'etp_left_Monitor', cv);
    [~, Etp_r]  = load_data_robust(folder, 'etp_right_Monitor', cv);

    [~, Xth_b]  = load_data_robust(folder, 'xth_bot_Monitor', cv);
    [~, Xth_t]  = load_data_robust(folder, 'xth_top_Monitor', cv);
    [~, Xth_l]  = load_data_robust(folder, 'xth_left_Monitor', cv);
    [~, Xth_r]  = load_data_robust(folder, 'xth_right_Monitor', cv);

    dek_dt = gradient(ek, t);
    dxv_dt = gradient(mass, t);

    xm  = Etp_b + Etp_t + Etp_l + Etp_r;
    xv  = dek_dt + dxv_dt;
    xth = Xth_b + Xth_t + Xth_l + Xth_r;
    xk  = xm + xv;
    tot = aphi + anab + xm + xth + xv;

    F.t    = t;
    F.aphi = aphi / power_ref;
    F.anab = anab / power_ref;
    F.xm   = xm   / power_ref;
    F.xth  = xth  / power_ref;
    F.xv   = xv   / power_ref;
    F.xk   = xk   / power_ref;
    F.tot  = tot  / power_ref;
end

% -- Draws one full Exergy decomposition (all terms) into a given axes,
%    used to populate each CV tab in the two tab-group figures. Styling
%    (colors, line width, latex labels/legend) matches the Exergy Balance
%    plot (Method 2) in drag_compare_STAR_4_fix_step_ALEX.m exactly. The
%    x-axis is normalized by the shedding period (t/T), starting at
%    t_eval_start, showing exactly num_cycles periods -- same logic as
%    time_norm/eval_mask in drag_compare_STAR_4_fix_step_ALEX.m. --
function plot_decomp_tab(ax, data, term_labels_latex, term_colors, term_linestyles, W, cv, ...
                          solver_name, t_eval_start, T_shed, num_cycles)
    fields = {'aphi', 'anab', 'xm', 'xth', 'xv', 'xk', 'tot'};

    time_norm = (data.t - t_eval_start) / T_shed;
    mask = (time_norm >= 0) & (time_norm <= num_cycles);
    t_plot = time_norm(mask);

    hold(ax, 'on');
    for k = 1:numel(fields)
        y_plot = data.(fields{k})(mask);
        plot(ax, t_plot, y_plot, 'Color', term_colors{k}, 'LineStyle', term_linestyles{k}, ...
             'LineWidth', W.line, 'DisplayName', term_labels_latex{k});
    end

    title(ax, sprintf('%s -- CV%d Exergy Balance', solver_name, cv), 'Interpreter', 'latex');
    xlabel(ax, '$t/T$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(ax, 'Exergy Coefficient', 'FontSize', 12);
    xlim(ax, [0, num_cycles]);
    grid(ax, 'on'); ax.GridAlpha = 0.3; ax.XTick = 0:num_cycles;
    legend(ax, 'Location', 'southoutside', 'Interpreter', 'latex', 'NumColumns', 4, 'FontSize', 14, 'Box', 'off');
end

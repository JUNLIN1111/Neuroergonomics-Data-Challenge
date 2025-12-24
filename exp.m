%% =========================================================
%  MRes Neuroergonomics – Data Challenge (Main Task)
%  Full analysis script (Segregational + Integrational)
%  Author: <Your Name>
% =========================================================

clear; clc; close all;
%% ===================== Channel definitions =====================

% Left / Right hemisphere (example mapping)
left_channels  = [1 2 3 4 5 6 7 8 9 10 11 12];
right_channels = [13 14 15 16 17 18 19 20 21 22 23 24];

% Orbital (inferior PFC) vs Lateral PFC
orbital_channels = [5 6 11 12 17 18 23 24];
lateral_channels = setdiff(1:24, orbital_channels);

% Sanity check
fprintf('Left channels: %s\n', mat2str(left_channels));
fprintf('Right channels: %s\n', mat2str(right_channels));
fprintf('Orbital channels: %s\n', mat2str(orbital_channels));
fprintf('Lateral channels: %s\n', mat2str(lateral_channels));

%% ===================== 1. Load database ==================
% !!! Modify to your own file path !!!
% Note: If the file does not exist, check the path or create test data first
% DB = readtable('D:\Machine_learing\data\MainTask\SPNDataChallenge0001.csv');

% [For testing] Create simulated dataset (Uncomment below to run without real data)
% ==========================================================
n_rows = 1000;
DB = table(...
    randi([1,20],n_rows,1), ...                  % SUBJECT
    randi([1,3],n_rows,1), ...                  % SESSION (1=Novice,2=Registrar,3=Consultant)
    repmat({'fNIRS'},n_rows,1), ...             % DATASOURCE
    repmat({''},n_rows,1), ...                  % STRUCTUREDDATA
    randi([1,24],n_rows,1), ...                 % CHANNEL
    randi([1,4],n_rows,1), ...                  % SIGNAL (1=OXY-HB)
    randi([1,5],n_rows,1), ...                  % STIMULUS
    randi([1,10],n_rows,1), ...                 % BLOCK
    randn(n_rows,1)*2 + 5, ...                  % MEAN_BASELINE
    rand(n_rows,1)*1 + 0.5, ...                 % STD_BASELINE
    randn(n_rows,1)*2 + 7, ...                  % MEAN_TASK
    rand(n_rows,1)*1 + 0.5, ...                 % STD_TASK
    randn(n_rows,1)*2.5, ...                    % TASK_MINUS_BASELINE
    randn(n_rows,1)*10 + 50, ...                % AREA_UNDER_CURVE_TASK
    rand(n_rows,1)*5 + 1, ...                   % TIME_TO_PEAK
    rand(n_rows,1)*5 + 1 ...                    % TIME_TO_NADIR
);
DB.Properties.VariableNames = { ...
    'SUBJECT', 'SESSION', 'DATASOURCE', 'STRUCTUREDDATA', ...
    'CHANNEL', 'SIGNAL', 'STIMULUS', 'BLOCK', ...
    'MEAN_BASELINE', 'STD_BASELINE', 'MEAN_TASK', 'STD_TASK', ...
    'TASK_MINUS_BASELINE', 'AREA_UNDER_CURVE_TASK', ...
    'TIME_TO_PEAK', 'TIME_TO_NADIR'};
% ==========================================================

disp('Database loaded successfully.');
disp(size(DB));
disp(DB.Properties.VariableNames);

%% ===================== 2. Select HbO ======================
% Signal identifiers:
% 1 = OXY-HB, 2 = DEOXY-HB, 4 = TOTAL-HB
HB_OXY = 1;
DB = DB(DB.SIGNAL == HB_OXY, :);

%% ===================== 3. Extract group ===================
% Expertise level encoded in SESSION:
% 1 = Novice, 2 = Registrar, 3 = Consultant
Group = strings(height(DB),1);
Group(DB.SESSION == 1) = "Novice";
Group(DB.SESSION == 2) = "Registrar";
Group(DB.SESSION == 3) = "Consultant";
DB.Group = categorical(Group);

disp('Groups detected:');
disp(categories(DB.Group));

%% ===================== 4. Define ROIs =====================
% Nominal ROI grouping (as allowed in coursework)
ROI.left    = [1 2 3 4 5 6];
ROI.right   = [7 8 9 10 11 12];
ROI.orbital = [13 14 15 16];
ROI.lateral = 1:12;

%% =========================================================
%% Hypothesis 1 – Segregational analysis
% Lateralization: Left vs Right
%% =========================================================
fprintf('\n=== Detailed stats: Lateralization ===\n');
groups = categories(DB.Group);

for g = 1:numel(groups)
    grp = groups{g};
    idx = DB.Group == grp & DB.SIGNAL == 1; % HbO

    left_idx  = idx & ismember(DB.CHANNEL, left_channels);
    right_idx = idx & ismember(DB.CHANNEL, right_channels);

    left_vals  = DB.TASK_MINUS_BASELINE(left_idx);
    right_vals = DB.TASK_MINUS_BASELINE(right_idx);

    nL = numel(left_vals);
    nR = numel(right_vals);

    mL = mean(left_vals);  sL = std(left_vals);
    mR = mean(right_vals); sR = std(right_vals);

    % Cohen's d (Handle cases with zero sample size)
    if nL>1 && nR>1
        sp = sqrt(((nL-1)*sL^2 + (nR-1)*sR^2)/(nL+nR-2));
        d  = (mL - mR)/sp;
    else
        d = NaN;
    end

    fprintf('%s\n', grp);
    fprintf('  Left : n=%d, mean=%.4f, std=%.4f\n', nL, mL, sL);
    fprintf('  Right: n=%d, mean=%.4f, std=%.4f\n', nR, mR, sR);
    fprintf('  Cohen''s d = %.3f\n\n', d);
end

%% =========================================================
%% Hypothesis 2 – Novices: Orbital > Lateral
%% =========================================================
fprintf('\n=== Detailed stats: Novice Orbital vs Lateral ===\n');
idx = DB.Group == "Novice" & DB.SIGNAL == 1;

orb_idx = idx & ismember(DB.CHANNEL, orbital_channels);
lat_idx = idx & ismember(DB.CHANNEL, lateral_channels);

orb_vals = DB.TASK_MINUS_BASELINE(orb_idx);
lat_vals = DB.TASK_MINUS_BASELINE(lat_idx);

nO = numel(orb_vals);
nL = numel(lat_vals);

mO = mean(orb_vals); sO = std(orb_vals);
mL = mean(lat_vals); sL = std(lat_vals);

% Cohen's d (Handle cases with zero sample size)
if nO>1 && nL>1
    sp = sqrt(((nO-1)*sO^2 + (nL-1)*sL^2)/(nO+nL-2));
    d  = (mO - mL)/sp;
else
    d = NaN;
end

fprintf('Orbital : n=%d, mean=%.4f, std=%.4f\n', nO, mO, sO);
fprintf('Lateral : n=%d, mean=%.4f, std=%.4f\n', nL, mL, sL);
fprintf('Cohen''s d = %.3f\n', d);

%% =========================================================
%% Hypothesis 3 – Integrational analysis
% "Prefrontal response is sensitive to surgical expertise"
%% =========================================================
fprintf('\n=== Hypothesis 3: Expertise effect (ANOVA) ===\n');
Y = DB.TASK_MINUS_BASELINE;
G = DB.Group;

[p,tbl,stats] = anova1(Y, G, 'off');
fprintf('One-way ANOVA: p = %.4f\n', p);

if p < 0.05 && ~isnan(p)
    fprintf('Post-hoc multiple comparisons:\n');
    multcompare(stats);
end

fprintf('\n=== Group-wise descriptive stats (Expertise) ===\n');
for g = 1:numel(groups)
    grp = groups{g};
    vals = DB.TASK_MINUS_BASELINE(DB.Group==grp & DB.SIGNAL==1);
    fprintf('%s: n=%d, mean=%.4f, std=%.4f\n', ...
        grp, numel(vals), mean(vals), std(vals));
end

%% ===================== 5. Figures =========================
% Vibrant color scheme (Ensure RGB values are between 0 and 1)
novice_color    = [0.9290  0.6940  0.1250];   % Orange yellow - Novice
registrar_color = [0.3010  0.7450  0.9330];   % Sky blue - Registrar
consultant_color= [0.6350  0.0780  0.1840];   % Dark red - Consultant
orb_color       = [0.2  0.8  0.4];            % Bright green - Orbitofrontal
lat_color       = [0.6  0.2  0.8];            % Purple - Lateral

%% --- Figure 1: Expertise effect ---
figure('Position', [100, 100, 900, 600]);
hold on;
groups = categories(DB.Group);
colors = [novice_color; registrar_color; consultant_color];

for i = 1:numel(groups)
    grp = groups{i};
    idx = DB.Group == grp;
    % Assign x-axis position i for each group
    boxchart(i*ones(sum(idx),1), DB.TASK_MINUS_BASELINE(idx), 'BoxFaceColor', colors(i,:), ...
             'MarkerStyle', 'o', 'MarkerSize', 5, 'LineWidth', 2.5);
end
hold off;

set(gca, 'FontSize', 16);
xlabel('Surgical Expertise Level', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('HbO Task − Baseline (μM)', 'FontSize', 20, 'FontWeight', 'bold');
title('Prefrontal HbO Activation by Expertise Level', ...
    'FontSize', 24, 'FontWeight', 'bold');
grid on; box on; set(gca, 'LineWidth', 1.5);
xticks(1:numel(groups));
xticklabels(groups);
xlim([0.5 numel(groups)+0.5]);

%% --- Figure 2: Novice orbital vs lateral ---
figure('Position', [100, 100, 800, 600]);
hold on;

regions = {'Orbitofrontal', 'Lateral'};
colors = [orb_color; lat_color];

% Assign x-axis positions 1 and 2 for each region
for i = 1:numel(regions)
    if i == 1
        vals = orb_vals;
    else
        vals = lat_vals;
    end
    x = i * ones(numel(vals),1);
    boxchart(x, vals, 'BoxFaceColor', colors(i,:), ...
             'MarkerStyle', 'o', 'MarkerSize', 5, 'LineWidth', 2.5);
end

hold off;

% Bold median line (white for better visibility, default Matlab tag)
hMed = findobj(gca, 'Tag', 'Median');
set(hMed, 'LineWidth', 3, 'Color', [1 1 1]);

% Set outliers to gray
hOut = findobj(gca, 'Tag', 'Outlier');
if ~isempty(hOut)
    set(hOut, 'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerFaceColor', [0.6 0.6 0.6]);
end

set(gca, 'FontSize', 16);
xlabel('Prefrontal Region', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('HbO Task − Baseline (μM)', 'FontSize', 20, 'FontWeight', 'bold');
title('Novice Group: Orbitofrontal vs. Lateral Activation', ...
    'FontSize', 24, 'FontWeight', 'bold');
grid on; box on; set(gca, 'LineWidth', 1.5);

xticks(1:numel(regions));
xticklabels(regions);
xlim([0.5 numel(regions)+0.5]);

legend(regions, 'FontSize', 16, 'Location', 'best');

%% --- Figure 3: Area Under Curve by expertise ---
figure('Position', [100, 100, 900, 600]);
hold on;
groups = categories(DB.Group);
colors = [novice_color; registrar_color; consultant_color];

for i = 1:numel(groups)
    grp = groups{i};
    idx = DB.Group == grp;
    x = i * ones(sum(idx),1);
    boxchart(x, DB.AREA_UNDER_CURVE_TASK(idx), 'BoxFaceColor', colors(i,:), ...
             'MarkerStyle', 'o', 'MarkerSize', 5, 'LineWidth', 2.5);
end
hold off;

set(gca, 'FontSize', 16);
xlabel('Surgical Expertise Level', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Area Under Curve (Task)', 'FontSize', 20, 'FontWeight', 'bold');
title('HbO Response AUC by Surgical Expertise', ...
    'FontSize', 24, 'FontWeight', 'bold');
grid on; box on; set(gca, 'LineWidth', 1.5);
xticks(1:numel(groups));
xticklabels(groups);
xlim([0.5 numel(groups)+0.5]);

%% --- Figure 4: Time to Peak by expertise ---
figure('Position', [100, 100, 900, 600]);
hold on;
for i = 1:numel(groups)
    grp = groups{i};
    idx = DB.Group == grp;
    x = i * ones(sum(idx),1);
    boxchart(x, DB.TIME_TO_PEAK(idx), 'BoxFaceColor', colors(i,:), ...
             'MarkerStyle', 'o', 'MarkerSize', 5, 'LineWidth', 2.5);
end
hold off;
set(gca, 'FontSize', 16);
xlabel('Surgical Expertise Level', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Time to Peak (s)', 'FontSize', 20, 'FontWeight', 'bold');
title('HbO Time to Peak by Surgical Expertise', 'FontSize', 24, 'FontWeight', 'bold');
grid on; box on; set(gca, 'LineWidth', 1.5);
xticks(1:numel(groups));
xticklabels(groups);
xlim([0.5 numel(groups)+0.5]);

%% --- Figure 5: Time to Nadir by expertise ---
figure('Position', [100, 100, 900, 600]);
hold on;
for i = 1:numel(groups)
    grp = groups{i};
    idx = DB.Group == grp;
    x = i * ones(sum(idx),1);
    boxchart(x, DB.TIME_TO_NADIR(idx), 'BoxFaceColor', colors(i,:), ...
             'MarkerStyle', 'o', 'MarkerSize', 5, 'LineWidth', 2.5);
end
hold off;
set(gca, 'FontSize', 16);
xlabel('Surgical Expertise Level', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Time to Nadir (s)', 'FontSize', 20, 'FontWeight', 'bold');
title('HbO Time to Nadir by Surgical Expertise', 'FontSize', 24, 'FontWeight', 'bold');
grid on; box on; set(gca, 'LineWidth', 1.5);
xticks(1:numel(groups));
xticklabels(groups);
xlim([0.5 numel(groups)+0.5]);

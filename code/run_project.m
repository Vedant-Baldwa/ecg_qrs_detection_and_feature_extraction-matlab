% run_project.m
clc; clear; close all;
project_root = pwd;
data_dir = fullfile(project_root,'data');
code_dir = fullfile(project_root,'code');
results_dir = fullfile(project_root,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

fprintf('Running ECG project...\n');

% --- Load CSV (mitbih_train.csv expected) ---
csvfile = fullfile(data_dir,'mitbih_train.csv');
if ~exist(csvfile,'file')
    error('Put mitbih_train.csv into the data/ folder and re-run.');
end
M = readmatrix(csvfile);
[nrows,ncols] = size(M);
siglen = ncols - 1;
snips = double(M(:,1:siglen));
labels = M(:,end);

% --- Detection (calls detect_r_peaks) ---
params.fs = 360;            % nominal sampling (used for filter design)
params.bp_low = 5; params.bp_high = 25; params.bp_order = 2;
params.prom_scale = 0.25;   % per-snippet prominence fraction (0.1-0.5)
params.min_dist_ms = 120;   % min distance between peaks (ms)
params.center_window_ms = 80; % prefer peaks within +/- this around center
% run detection
[r_locs, r_vals, peak_found] = detect_r_peaks(snips, params);

% Save detection
save(fullfile(results_dir,'final_detection.mat'),'r_locs','r_vals','peak_found','labels','params','-v7.3');

% --- Analysis & plots (calls analyze_and_plot) ---
analyze_and_plot(snips, labels, r_locs, results_dir, params);

fprintf('Done. Results & figures are in %s\n', results_dir);

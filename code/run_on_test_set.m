%% run_on_test_set.m
% Run the minimal detector on mitbih_test.csv and compute central-window accuracy.

project_root = pwd;
data_dir = fullfile(project_root,'data');
results_dir = fullfile(project_root,'results');
test_out = fullfile(results_dir,'test');
if ~exist(test_out,'dir'), mkdir(test_out); end

% 1) Load test CSV
test_csv = fullfile(data_dir,'mitbih_test.csv');
if ~exist(test_csv,'file')
    error('mitbih_test.csv not found in data/. Place it there and re-run.');
end
Mtest = readmatrix(test_csv);
[nrows_test, ncols_test] = size(Mtest);
siglen_test = ncols_test - 1;
snips_test = double(Mtest(:,1:siglen_test));
labels_test = Mtest(:,end);

% 2) Load params from your training run (saved in results_final/final_detection.mat)
train_mat = fullfile(results_dir,'final_detection.mat');
if exist(train_mat,'file')
    S = load(train_mat);
    if isfield(S,'params')
        params = S.params;
        fprintf('Loaded params from %s\n', train_mat);
    else
        % fallback default params (should match your run_project.m)
        params.fs = 360; params.bp_low = 5; params.bp_high = 25; params.bp_order = 2;
        params.prom_scale = 0.25; params.min_dist_ms = 120; params.center_window_ms = 80;
        fprintf('Params not found in %s — using fallback defaults.\n', train_mat);
    end
else
    error('Training results not found: %s. Run run_project.m on train set first.', train_mat);
end

% 3) Run detection (uses your detect_r_peaks.m from code/)
addpath(fullfile(project_root,'code')); % ensure function on path
[r_locs_t, r_vals_t, peak_found_t] = detect_r_peaks(snips_test, params);

% 4) Save raw detection
save(fullfile(test_out,'test_detection.mat'),'r_locs_t','r_vals_t','peak_found_t','labels_test','params','-v7.3');

% 5) Compute central-window accuracy (±20 samples) and per-label stats
centerIdx = ceil(siglen_test/2);
win = 20;
overall = mean(abs(r_locs_t - centerIdx) <= win) * 100;

u = unique(labels_test);
T = table('Size',[numel(u) 2],'VariableTypes',{'double','double'},...
    'VariableNames',{'Label','CentralAcc'});
for k = 1:numel(u)
    idx = labels_test==u(k);
    T.Label(k) = u(k);
    T.CentralAcc(k) = mean(abs(r_locs_t(idx)-centerIdx) <= win) * 100;
end

% 6) Save & plot simple histograms for comparison with training set
writetable(T, fullfile(test_out,'central_window_accuracy_test.csv'));
fid = fopen(fullfile(test_out,'summary.txt'),'w');
fprintf(fid,'CentralWindow = %d\nOverallCentralAcc_percent = %.4f\n\n', win, overall);
fclose(fid);

fprintf('\nTest set results (central window = ±%d samples):\n', win);
fprintf('OverallCentralAcc_percent = %.4f\n', overall);
disp(T);

% Histogram of R locations
figure('Visible','off'); histogram(r_locs_t, min(40,siglen_test)); title('Test: R location histogram');
xlabel('R sample index'); ylabel('Count'); saveas(gcf, fullfile(test_out,'r_loc_hist_test.png'));

% Per-label histogram (save)
figure('Visible','off','Position',[100 100 1200 800]);
numLabels = numel(u);
for k = 1:numLabels
    lab = u(k);
    idxs = find(labels_test==lab);
    subplot(ceil(numLabels/2),2,k);
    histogram(r_locs_t(idxs), min(30,siglen_test));
    title(sprintf('Label %g (n=%d)', lab, numel(idxs)));
    xlabel('R sample'); ylabel('Count');
end
saveas(gcf, fullfile(test_out,'r_loc_hist_by_label_test.png'));

% Copy one-line summary to console for easy paste into report
fprintf('\nSaved test results to %s\n', test_out);

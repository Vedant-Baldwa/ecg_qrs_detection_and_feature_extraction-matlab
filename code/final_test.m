load('results/final_detection.mat','r_locs','labels','params');
M = readmatrix('data/mitbih_train.csv');
L = size(M,2)-1;
centerIdx = ceil(L/2);
win = 20;  % try 20 samples (~55 ms if fs=360)

overall = mean(abs(r_locs - centerIdx) <= win) * 100;

u = unique(labels);
T = table('Size',[numel(u) 2],'VariableTypes',{'double','double'},...
    'VariableNames',{'Label','CentralAcc'});
for k = 1:numel(u)
    idx = labels==u(k);
    T.Label(k) = u(k);
    T.CentralAcc(k) = mean(abs(r_locs(idx)-centerIdx) <= win) * 100;
end

disp(struct('CentralWindow',win,'OverallCentralAcc_percent',overall));
disp(T)
writetable(T, 'results/central_window_accuracy.csv');
function analyze_and_plot(snips, labels, r_locs, results_dir, params)
% analyze_and_plot - creates required plots, templates and prints stats

if ~exist(results_dir,'dir'), mkdir(results_dir); end
[N,L] = size(snips);
centerIdx = ceil(L/2);

% Basic stats
found_pct = sum(~isnan(r_locs))/N*100;
fprintf('Detection: %d / %d (%.2f%%) snippets\n', sum(~isnan(r_locs)), N, found_pct);

% R-location histogram overall
figure('Visible','off'); histogram(r_locs, min(40,L));
xlabel('R location (sample index)'); ylabel('Count'); title('R location histogram (all snippets)');
saveas(gcf, fullfile(results_dir,'r_loc_hist_all.png'));

% Per-label histograms and aligned mean templates
uniq = unique(labels);
figure('Visible','off','Position',[100 100 1200 800]);
numLabels = numel(uniq);
for k = 1:numLabels
    lab = uniq(k);
    idxs = find(labels==lab);
    locs = r_locs(idxs);
    subplot(ceil(numLabels/2),2,k);
    histogram(locs, min(30,L));
    title(sprintf('Label %g (n=%d)',lab,numel(idxs)));
    xlabel('R sample'); ylabel('Count');
end
saveas(gcf, fullfile(results_dir,'r_loc_hist_by_label.png'));

% Aligned mean templates per label
figure('Visible','off','Position',[100 100 1200 800]);
for k = 1:numLabels
    lab = uniq(k);
    idxs = find(labels==lab);
    mat = zeros(numel(idxs), L);
    cnt = 0;
    for j = 1:numel(idxs)
        i = idxs(j);
        r = r_locs(i);
        if isnan(r), continue; end
        row = snips(i,:)';
        % shift so that r -> center
        shift = centerIdx - round(r);
        out = zeros(L,1);
        inIdx = (1:L)'; outIdx = inIdx + shift;
        valid = outIdx>=1 & outIdx<=L;
        out(outIdx(valid)) = row(inIdx(valid));
        cnt = cnt + 1;
        mat(cnt,:) = out';
    end
    if cnt==0, mu = zeros(1,L); sd = zeros(1,L); else
        mat = mat(1:cnt,:);
        % normalize rows by max abs to compare shapes
        for rr=1:cnt
            row = mat(rr,:);
            m = max(abs(row));
            if m>0, mat(rr,:) = row / m; end
        end
        mu = mean(mat,1); sd = std(mat,0,1);
    end
    subplot(ceil(numLabels/2),2,k); hold on;
    x = (1:L) - centerIdx;
    fill([x fliplr(x)], [mu+sd fliplr(mu-sd)], [0.9 0.9 0.9], 'EdgeColor','none');
    plot(x, mu, 'b', 'LineWidth', 1.5); xline(0,'r--');
    title(sprintf('Label %g mean template (n=%d)', lab, cnt));
    xlabel('Samples relative to R (center=0)'); ylabel('Normalized amplitude');
    hold off;
end
saveas(gcf, fullfile(results_dir,'templates_by_label.png'));

% Save small CSV of summary stats
summary_table = table();
for k = 1:numLabels
    lab = uniq(k);
    idxs = find(labels==lab);
    locs = r_locs(idxs);
    locs = locs(~isnan(locs));
    summary_table = [summary_table; ...
        table(lab, numel(idxs), numel(locs), mean(locs), std(locs), 'VariableNames',{'Label','Count','Found','MeanLoc','StdLoc'})];
end
writetable(summary_table, fullfile(results_dir,'summary_stats_by_label.csv'));

% Save a few example snippets (first 12) for demo slide
demoDir = fullfile(results_dir,'demo_examples');
if ~exist(demoDir,'dir'), mkdir(demoDir); end
nDemo = min(12, N);
for i = 1:nDemo
    hf = figure('Visible','off'); plot(snips(i,:)); hold on; plot(r_locs(i), snips(i,r_locs(i)),'ro'); title(sprintf('Row %d label=%g', i, labels(i)));
    saveas(hf, fullfile(demoDir, sprintf('example_%02d.png', i)));
    close(hf);
end

fprintf('Analysis figures & summary saved to %s\n', results_dir);
end

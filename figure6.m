
plot_figure6()



function plot_figure6(options)
    arguments
        options.bandRange = [40 100];
        options.isGroupRecs = true;
    end
    nrows = 5;
    ncols = 4;
    fig = figure("Name", "ephys_figure");
    set(fig, 'Units','centimeters', 'Position',[2 2 12 17]); % 8.5 cm wide (single column)
    set(groot,'DefaultAxesFontName','Helvetica');
    set(groot,'DefaultTextFontName','Helvetica');
    set(groot, ...
        'defaultAxesFontSize', 8, ...                    % tick labels
        'defaultAxesLabelFontSizeMultiplier', 10/8, ...  % x/y labels
        'defaultAxesTitleFontSizeMultiplier', 12/8);     % titles
    tdl = tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','compact','TileIndexing','columnmajor');
    
    fprintf('\nStart plotting figure 6...\n')
    pb = predictionBreak;
    stats_csv = fullfile(pb.outputDir, 'figure6_stats.csv');
    if isfile(stats_csv)
        delete(stats_csv);
    end
    data = {
        {'PV106', 'PV157', 'PV153', 'PV126'}, 'escape_time', 'circle';
        % {'PV157', 'PV106', 'PV153', 'PV162'}, 'escape_time', 'low_horizontal_noise';
        {'PV106', 'PV143', 'PV153'}, 'flip_time', 'circle'
    };
    lfp_examples = {
        'PV106', 'Hunter17', 12, [];
        % 'PV157', 'Hunter72', 14;
        % 'PV106', 'Hunter17', 15;
        'PV143', 'Hunter16', 4, [-5 1];
    };
    spk = {
        'PV106', 'Hunter17', 12, [4 7 15 17];
        'PV157', '', 0, [35 43 51 54 87];
        'PV153', '', 0, [12 13 16 17 18 19 20 24 28 43 45 59 60];
        'PV126', 'Hunter45', 34, [10,11,12,22,45,50,51];
        % 'PV157', 'Hunter81', 25, [38 39 67 73 74 81 83];
        % 'PV106', '', 0, [68 69 73 77];
        % 'PV153', '', 0, [5 8 25 26 27 41 62];
        % 'PV162', '', 0, 48;
        'PV106', 'Hunter17', [12 15 18 29 30 33], [17 32];
        'PV143', '', 0, [20 21];
        'PV153', '', 0, [10 16 60];
    };
    axes = gobjects(nrows,ncols);
    si = 1;
    [animal_col, power_col, label_col, metric_col] = deal({}, [], {}, {});
    for i=1:height(data)
        % raw data example + gamma band 
        axes(1,i) = nexttile([1 2]);
        pb.plot_event_segments(lfp_examples{i,1}, lfp_examples{i,2}, data{i,2}, lfp_examples{i,3}, ax=true, time_before=2, time_after=2,...
                               bandRange=options.bandRange, isPutTrialLabels=false, is_scalebar=i==1, is_cache=true);
        % lfp spectrogram example
        axes(2,i) = nexttile([1 2]);
        [S, recNamesIdx] = pb.plot_all_event_spectrogram(lfp_examples{i,1}, data{i,3}, data{i,2}, ax=true, ...
                                          time_before=1, time_after=1, all_together=true, min_trials=1, ...
                                          is_engaged=true, freqRange=[1 150], is_cache=true, clim=lfp_examples{i,4});
        title('');
        if i < ncols
            delete(colorbar(gca));
        else
            cb = colorbar(gca);
            cb.Location = 'eastoutside';
        end
        if (~options.isGroupRecs); recNamesIdx = nan; end

        % gamma band pre vs. post swarmplot
        axes(3,i) = nexttile;
        [pVal, stats] = pb.plot_band_statistics(S, bandRange=options.bandRange, recNamesIdx=recNamesIdx, isPlotPvalue=false);
        fprintf('%s LFP-swarmplot paired ttest: t(%d) = %.1f, p = %.2e\n', lfp_examples{i,1}, stats.df, stats.tstat, pVal);

        % spikes raster plot example
        axes(4,i) = nexttile([1 2]);
        pb.plot_spikes_raster_plot(spk{si,1}, spk{si,2}, data{i,2}, cluster_id=spk{si,3}, ax=true);

        % spike rate single animal
        axes(5,i) = nexttile;
        strike_recs = arrayfun(@(x) sprintf('Hunter%d', x), spk{si,4}, 'UniformOutput', false);
        [avg_rate, t, ~] = pb.plot_avg_spikes_rate(spk{si,1}, strike_recs, data{i,2}, sec_before=1, sec_after=1, binw=50, ax=true, ...
                                                    smooth_kernel=11, is_cache=true);

        % lfp gamma comparisons for all animals
        animals = data{i,1};
        axes(6,i) = nexttile;
        for j=1:numel(animals)
            if j > 1  % no need to run again
                [S, recNamesIdx] = pb.plot_all_event_spectrogram(animals{j}, data{i,3}, data{i,2}, isPlot=false, ...
                                              time_before=1, time_after=1, all_together=true, min_trials=1, ...
                                              is_engaged=true, freqRange=[1 150], is_cache=true);
                if (~options.isGroupRecs || i>1); recNamesIdx = nan; end
            end
            [pre, post] = pb.get_band_means(S, bandRange=options.bandRange, recNamesIdx=recNamesIdx);
            plot_greater_than_zero_stat(post, pre, j, animals, 2.3)
            power_diff = post - pre;
            append_figure6_stats(stats_csv, animals{j}, 'lfp', data{i,2}, power_diff);
            animal_col = [animal_col; repmat({animals{j}}, numel(power_diff), 1)];
            label_col = [label_col; repmat({data{i,2}}, numel(power_diff), 1)];
            metric_col = [metric_col; repmat({'lfp'}, numel(power_diff), 1)];
            power_col  = [power_col; power_diff(:)];
        end
        set(gca,'XTick',1:numel(animals),'XTickLabel',animals);
        
        % spikes rate comparisons all animals
        axes(7,i) = nexttile;
        for j=1:numel(animals)
            if isempty(spk{si,4})
                si = si + 1;
                continue
            end
            strike_recs = arrayfun(@(x) sprintf('Hunter%d', x), spk{si,4}, 'UniformOutput', false);
            [avg_rate, t, ~] = pb.get_spike_rates(animals{j}, strike_recs, data{i,2}, is_engaged=true, ...
                                                    sec_before=1, sec_after=1, binw=50, ...
                                                    is_tuned=true, alpha=0.05, is_cache=true);
            if isempty(avg_rate)
                fprintf('No spike-rate data for %s, %s; skipping statistics panel.\n', animals{j}, data{i,2});
                si = si + 1;
                continue
            end
            pre = mean(avg_rate(:,t<=0),2);
            post = mean(avg_rate(:,t>0),2);
            plot_greater_than_zero_stat(post, pre, j, animals, nan)
            rate_diff = post - pre;
            append_figure6_stats(stats_csv, animals{j}, 'spike_rate', data{i,2}, rate_diff);
            animal_col = [animal_col; repmat({animals{j}}, numel(rate_diff), 1)];
            label_col = [label_col; repmat({data{i,2}}, numel(rate_diff), 1)];
            metric_col = [metric_col; repmat({'spike_rate'}, numel(rate_diff), 1)];
            power_col  = [power_col; rate_diff(:)];
            si = si + 1;
        end
    end
    
    df = table(animal_col, metric_col, label_col, power_col, 'VariableNames', {'animal_id', 'metric', 'label', 'power_diff'});
    
    for r = [1 2 3 6]
        linkaxes(axes(r,1:end-1), 'y');
    end
    % for c = 2:(ncols-1)
    %     set(axes(1:4,c), 'YTickLabel', []);
    %     ylabel(axes(1:4, c), '')
    % end
    % set(findall(fig,'Type','line'),'LineWidth',0.75);
    set(findall(fig,'Type','scatter'),'SizeData',10); % adjust if needed
    exportgraphics(fig, sprintf('%s/fig6_script.pdf', pb.outputDir), 'ContentType','vector')
end


function plot_greater_than_zero_stat(post, pre, j, animals, ymax)
    dp = post-pre;
    [~,p,ci,stats] = ttest(post, pre, 'Tail','right');
    d = mean(dp)/std(dp,0);                % Cohen's d
    g = d*(1 - 3/(4*numel(dp) - 9));      % Hedges' g (small-sample corrected)
    % h_norm = lillietest(dp);
    ste = std(dp, 0, 'omitnan')/sqrt(sum(isfinite(dp)));
    fprintf('%s: t(%d) = %.1f, p = %.3e, CohenHedges: %.3f, %.1f±%.1f\n', animals{j}, stats.df, stats.tstat, p, g, mean(dp), ste);
    bar(j, mean(dp))
    hold on
    errorbar(j, mean(dp), ste, 'k')
    if ~isnan(ymax)
        ylim([0 ymax])
    else
        yl = ylim;
        ymax = yl(2);
    end
    text(j, mean(dp)+ste, p2stars(p),'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    box(gca,'off');
    % ylabel('Post-Pre');
end


function append_figure6_stats(csvFile, animalId, metric, label, values)
    values = values(:);
    if isempty(values)
        return
    end
    T = table(repmat({animalId}, numel(values), 1), ...
              repmat({metric}, numel(values), 1), ...
              repmat({label}, numel(values), 1), ...
              values, ...
              'VariableNames', {'animal_id', 'metric', 'label', 'power_diff'});
    if isfile(csvFile)
        writetable(T, csvFile, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(T, csvFile);
    end
end

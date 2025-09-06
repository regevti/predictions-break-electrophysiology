
plot_figure6()



function plot_figure6(options)
    arguments
        options.bandRange = [40 100];
        options.isGroupRecs = false;
    end

    pb = predictionBreak;
    figure("Name", "ephys_figure");
    data = {
        'PV106', 'escape_time', 'circle', 'Hunter17', 12, {'PV106', 'PV157', 'PV126'};
        'PV157', 'escape_time', 'low_horizontal_noise', 'Hunter72', 14, {'PV157', 'PV162'};
        'PV106', 'flip_time', 'circle', 'Hunter17', 15, {'PV106', 'PV143'}
        };
    for i=1:height(data)
        subplot(3, 4, 4*i-3)
        pb.plot_event_segments(data{i,1}, data{i,4}, data{i,2}, data{i,5}, ax=true, ...
                               bandRange=options.bandRange, isPutTrialLabels=false);
        subplot(3, 4, 4*i-2)
        [S, recNamesIdx] = pb.plot_all_event_spectrogram(data{i,1}, data{i,3}, data{i,2}, ax=true, ...
                                          time_before=1, time_after=1, all_together=true, ...
                                          is_engaged=true, freqRange=[1 150], is_cache=true);
        if (~options.isGroupRecs); recNamesIdx = nan; end
        subplot(3, 4, 4*i-1)
        pb.plot_band_statistics(S, bandRange=options.bandRange, recNamesIdx=recNamesIdx);
            
        animals = data{i,6};
        subplot(3, 4, 4*i)
        for j=1:numel(animals)
            if j > 1  % no need to run again
                [S, recNamesIdx] = pb.plot_all_event_spectrogram(animals{j}, data{i,3}, data{i,2}, isPlot=false, ...
                                              time_before=1, time_after=1, all_together=true, min_trials=1, ...
                                              is_engaged=true, freqRange=[1 150], is_cache=true);
                if (~options.isGroupRecs); recNamesIdx = nan; end
            end
            [pre, post] = pb.get_band_means(S, bandRange=options.bandRange, recNamesIdx=recNamesIdx);
            dp = post-pre;
            bar(j, mean(dp))
            hold on
            errorbar(j, mean(dp), std(dp, 0, 'omitnan')/sqrt(sum(isfinite(dp))), 'k')
            title('Post-Pre')
        end
        set(gca,'XTick',1:numel(animals),'XTickLabel',animals);
    end
end


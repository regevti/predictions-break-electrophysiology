classdef predictionBreak
    %predictionBreak Summary of this class goes here
    %   Detailed explanation goes here

    properties
        wa
        mainT
        frameGap
        outputDir
        cacheDir
    end

    methods
        function obj = predictionBreak()
            % constructor
            fclose('all');
            obj.wa = wakeAnalysis('~/PhD/Matlab/brainStatesWake_mac.xlsx');
            mainT_path = '~/PhD/Matlab/prediction_break_trials.csv';
            opts = detectImportOptions(mainT_path);
            opts = setvartype(opts,'flip_time','datetime');
            obj.mainT = readtable(mainT_path, opts);
            obj.frameGap = readtable('~/PhD/Matlab/frames_gap.csv');
            obj.outputDir = '~/PhD/Matlab/plots';
            obj.cacheDir = '~/PhD/Matlab/cache';
        end


        function [arena, T, channel] = load_rec(obj, animalId, recName)
            % if ~isempty(obj.wa.currentDataObj)
            %     obj.wa.currentDataObj.closeOpenFiles();
            % end
            obj.wa.setCurrentRecording(sprintf('Animal=%s,recNames=%s', animalId, recName));
            obj.wa.getFilters();
            arena = obj.wa.getArenaCSVs(1);
            % get default channel
            R = obj.wa.recTable;
            channel = R.defaulLFPCh(strcmp(R.Animal,animalId)&(strcmp(R.recNames,recName)));
            
            T = obj.mainT(strcmp(obj.mainT.animal_id, animalId) & strcmp(obj.mainT.rec_name, recName), :);
        end


        function res = load_event_signals(obj, animalId, movementType, eventName, options)
            arguments
                obj predictionBreak
                animalId char
                movementType char {mustBeMember(movementType, ['circle','low_horizontal_noise'])}
                eventName char {mustBeMember(eventName, ['enter_time','escape_time', 'exit_time', 'flip_time'])}
                options.channel double = nan
                options.time_before double = 2 % time in sec
                options.time_after double = 2 % time in sec
                options.min_trials double = nan
                options.is_engaged logical = true
                options.no_strikes logical = true
                options.is_cache logical = true
            end
            cacheFile = obj.get_cache_file_name('event_signals', options, {'is_cache'}, animalId, movementType, eventName);
            if options.is_cache && isfile(cacheFile)
                load(cacheFile, 'res');
                return;
            end

            rec_names = obj.get_relevant_rec_names(animalId, eventName, movementType, ...
                                               is_engaged=options.is_engaged, no_strikes=options.no_strikes, ...
                                               min_trial=options.min_trials);
            res = cell(4, numel(rec_names)+1);
            for i=1:numel(rec_names)
                if ~mod(i,20)
                    fclose('all')
                end
                try
                    fprintf('>> (%d/%d) %s,%s\n', i, numel(rec_names), animalId, rec_names{i})
                    [arena, ~, defaultCh] = obj.load_rec(animalId, rec_names{i});
                    T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, rec_names{i}, '', options.no_strikes);
                    channel = defaultCh.*(isnan(options.channel)) + options.channel.*(~isnan(options.channel));
                    if isnan(options.channel), channel = defaultCh; end
                    [tEvent, ~] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, rec_names{i});
                    [V, ~] = obj.wa.currentDataObj.getData(channel, tEvent-(options.time_before*1000), (options.time_before+options.time_after)*1000);
                    V = squeeze(V);
                    if size(V, 2) == 1 % squeeze flip order if only one segment
                        V = V';
                    end
                    
                    fs = obj.wa.currentDataObj.samplingFrequency(channel);                    
                    res{1, i} = rec_names{i};
                    res{2, i} = V;
                    res{3, i} = fs;
                    res{4, i} = channel;
                catch ME
                    fprintf('ERROR: %s,%s; %s\n\n', animalId, rec_names{i}, ME.message);
                end
            end
            keep = ~cellfun('isempty', res(2,:));
            res  = res(:, keep);
            save(cacheFile, 'res');
        end

        
        function [S, recNamesIdx] = plot_all_event_spectrogram(obj, animalId, movementType, eventName, options)
            arguments
                obj predictionBreak
                animalId char
                movementType char {mustBeMember(movementType, ['circle','low_horizontal_noise'])}
                eventName char {mustBeMember(eventName, ['enter_time','escape_time', 'exit_time', 'flip_time'])}
                options.channel double = nan
                options.time_before double = 2 % time in sec
                options.time_after double = 2 % time in sec
                options.all_together logical = false
                options.freqRange double = [1 100]
                options.WindowSec double = 0.2
                options.StepSec double = 0.01
                options.min_trials double = nan
                options.is_engaged logical = true
                options.no_strikes logical = true
                options.is_cache logical = true
                options.isPlot logical = true
                options.ax logical = false % if true, figure is not created
                options.max_plots_per_figure double = 12;
                options.clim = []
            end
            
            res = obj.load_event_signals(animalId, movementType, eventName, channel=options.channel, time_before=options.time_before, ...
                                         time_after=options.time_after, min_trials=options.min_trials, is_engaged=options.is_engaged, ...
                                         no_strikes=options.no_strikes, is_cache=options.is_cache);
            heights = cellfun(@(x) size(x,1), res(2,:));
            recNamesIdx = repelem(1:numel(res(2,:)), heights);

            if ~options.ax && options.isPlot
                figure();
            end

            if ~options.all_together
                tlo = tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');
                currentPage = 1;
                title(tlo, sprintf('%s - page %d', animalId, currentPage));
                n_plotted = 0;
                nexttile;
                for i=1:width(res)
                    V = res{2, i};
                    title = sprintf('%s - %s - %s - %s\nn_trials = %d', animalId, res{1, i}, movementType, eventName, size(V, 1));
                    S = obj.run_spectrogram(V, mode([res{3, i}]), title=title, ...
                                           freqRange=options.freqRange, WindowSec=options.WindowSec, StepSec=options.StepSec, ...
                                           time_before=option.time_before, time_after=options.time_after, clim=options.clim);
                    n_plotted = n_plotted + 1;
                    if n_plotted >= options.max_plots_per_figure
                        figure();
                        tlo = tiledlayout('flow', 'TileSpacing','compact', 'Padding','compact');
                        currentPage = currentPage + 1;
                        title(tlo, sprintf('%s - page %d', animalId, currentPage));
                        n_plotted = 0;
                    end
                    nexttile;
                end

            else  % plot one panel with an average of all traces
                V = vertcat(res{2, :});
                title = sprintf('%s - %s - %s - %s\nn_trials = %d', animalId, 'all', movementType, eventName, size(V, 1));
                S = obj.run_spectrogram(V, mode([res{3, :}]), title=title, isPlot=options.isPlot, ...
                                           freqRange=options.freqRange, WindowSec=options.WindowSec, StepSec=options.StepSec, ...
                                           time_before=options.time_before, time_after=options.time_after, clim=options.clim);
            end
        end

        
        function S = run_spectrogram(obj, V, fs, options)
            arguments
                obj predictionBreak
                V double
                fs double
                options.title char = '';
                options.time_before double = 2 % time in sec
                options.time_after double = 2 % time in sec
                options.WindowSec double = 0.2
                options.StepSec double = 0.01
                options.freqRange double = [1 100]
                options.isPlot logical = true
                options.ax logical = true % if false create new figure
                options.clim = []
            end
            S = event_aligned_spectrogram(V, fs, ...
                                         'FreqRange', options.freqRange, ...
                                         'WindowSec', options.WindowSec, ...
                                         'StepSec',   options.StepSec, ...
                                         'PreWindow', [-options.time_before 0], ...
                                         'FullWindow',[-options.time_before options.time_after], ...
                                         'Title', options.title, ...
                                         'Plot', options.isPlot, ...
                                         'Clim', options.clim, ...
                                         'Ax', true);
        end

        
        function [pre, post] = get_band_means(obj, S, options)
            arguments
                obj predictionBreak
                S struct
                options.bandRange double = [30 60]
                options.recNamesIdx double = nan
            end
            band = squeeze(mean(S.P_event_dB((options.bandRange(1)<=S.f)&(S.f<=options.bandRange(2)),:, :), 1));
            if ~isnan(options.recNamesIdx)
                band = splitapply(@(x) mean(x, 2), band, options.recNamesIdx);
            end
            pre = mean(band(S.t<=0, :), 1);
            post = mean(band(S.t>0, :), 1);
        end


        function plot_band_statistics(obj, S, options)
            arguments
                obj predictionBreak
                S struct
                options.bandRange double = [30 60]
                options.recNamesIdx double = nan 
            end
            [pre, post] = obj.get_band_means(S, recNamesIdx=options.recNamesIdx, bandRange=options.bandRange);
            pre  = pre(:);
            post = post(:);

            [h,p,~,s] = ttest(pre, post);
            stats = s;

            x = [ones(size(pre)); 2*ones(size(post))];
            y = [pre; post];
            
            % Build color matrix BEFORE masking
            cPre  = [0.30 0.60 0.90];
            cPost = [0.90 0.40 0.40];
            C = [repmat(cPre,  numel(pre),  1);
                 repmat(cPost, numel(post), 1)];   % size must be numel(y) × 3
            
            % Mask out non-finite points from *both* y/x and colors
            m = isfinite(y) & isfinite(x);   % (x is finite anyway, but harmless)
            x = x(m);
            y = y(m);
            C = C(m, :);
            
            ax = gca; hold(ax,'on');
            swarmchart(ax, x, y, 36, C, 'filled');   % now rows(C) == numel(y)

            % means (optional visual)
            m1 = mean(pre);  m2 = mean(post);
            plot(ax, [1 2], [m1 m2], 'k-', 'LineWidth', 1.2);
            plot(ax, 1, m1, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
            plot(ax, 2, m2, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);

            % axes cosmetics
            xlim(ax, [0.5 2.5]);
            set(ax,'XTick',[1 2],'XTickLabel',{'pre','post'});
            ylabel(ax,'Signal');
            box(ax,'off');

            % significance bar + stars
            stars = p2stars(p);                      % '*', '**', '***' or 'n.s.'
            yMax  = max(y);
            yMin  = min(y);
            dy    = max(1e-12, range([yMin yMax]));
            yBar  = yMax + 0.06*dy;                  % bar height
            yTxt  = yBar + 0.03*dy;                  % text above bar

            plot(ax, [1 2], [yBar yBar], 'k-', 'LineWidth', 1.2);         % top bar
            plot(ax, [1 1], [yBar-0.02*dy yBar], 'k-', 'LineWidth', 1.2); % left tick
            plot(ax, [2 2], [yBar-0.02*dy yBar], 'k-', 'LineWidth', 1.2); % right tick
            txt = sprintf('p = %.3g  %s', p, stars);
            text(1.5, yTxt, txt, 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', 'FontWeight','bold');

            % expand ylim to fit annotation
            ylim(ax, [yMin-0.05*dy, yTxt+0.08*dy]);

            function s = p2stars(p)
                if p < 1e-3
                    s = '***';
                elseif p < 1e-2
                    s = '**';
                elseif p < 0.05
                    s = '*';
                else
                    s = 'n.s.';
                end
            end
        end
        

        function plot_event_segments(obj, animalId, recName, eventName, trialId, options)
            arguments
                obj predictionBreak
                animalId char
                recName char
                eventName char
                trialId double = nan
                options.is_engaged logical = false
                options.no_strikes = true
                options.time_before double = 2
                options.time_after double = 2
                options.bandPass double = nan
                options.channel = nan(1)
                options.ax logical = false
                options.bandRange double = nan
                options.isPutTrialLabels logical = true
            end
            [arena, ~, defaultCh] = obj.load_rec(animalId, recName);
            T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, recName);
            channel = defaultCh.*(isnan(options.channel)) + options.channel.*(~isnan(options.channel));
            if isnan(options.channel), channel = defaultCh; end
            if ~isnan(trialId)
                T = T(T.trial_id==trialId, :);
            end
            [tEvent, ~] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, recName);
            [V, t] = obj.wa.currentDataObj.getData(channel, tEvent-(options.time_before*1000), (options.time_before+options.time_after)*1000);
            fs = obj.wa.currentDataObj.samplingFrequency(channel);
            if ~isnan(options.bandPass)
                F = filterData(fs);
                F.lowPassCutoff = options.bandPass(1);
                F.highPassCutoff = options.bandPass(2);
                F = F.designBandPass;
                V = F.getFilteredData(V);
            end
            V = squeeze(V);
            if size(V, 2) == 1 % squeeze flip order if only one segment
                V = V';
            end

            if ~options.ax
                figure()
            end

            step = max(V,[],'all') - min(V,[],'all');
            for i=1:size(V, 1)
                hold on
                v = V(i, :);
                if ~isnan(options.bandRange)
                    vmin = min(v);  vmax = max(v);
                    v  = (v - vmin) ./ max(eps, vmax - vmin); 
                end
                plot(t/1000 - options.time_before, v+step*(i-1), 'k')
                if options.isPutTrialLabels
                    tagCol = sprintf('%s_tags', eventName);
                    trial_text = sprintf('trial%d', T.trial_id(i));
                    if ~options.is_engaged && strcmp(T.(tagCol){i}, 'not engaged')
                        trial_text = [trial_text ' (not_eng)'];
                    end
                    if ~options.no_strikes && strcmp(T.(tagCol){i}, 'strike')
                        trial_text = [trial_text ' (strike)'];
                    end
                    text(options.time_after+0.01, step*(i-1), trial_text)
                end
            end
            xline(0, 'r--')

            if ~isnan(options.bandRange)
                if size(V, 1) > 1; error('size of V must be 1'); end
                S = obj.run_spectrogram(V, fs, isPlot=false);
                band = mean(S.P_event_dB((options.bandRange(1)<=S.f)&(S.f<=options.bandRange(2)),:), 1);
                bmin = min(band); bmax = max(band);
                band  = (band - bmin) ./ max(eps, bmax - bmin);
                band = smoothdata(band, "sgolay", 17); 
                plot(S.t, (band/2)-0.5, 'Color', "#1171BE")
                scalebar(gca, 0, '', 5 / (bmax - bmin)/2, '5dB', 'sw');
                scalebar(gca, 0.1, '100ms', 100/(vmax - vmin), '100µV', 'nw');
                axis off
            end
            if options.isPutTrialLabels
                titleText = sprintf('%s, %s', animalId, recName);
                if ~isnan(options.bandPass)
                    titleText = [titleText sprintf('\nFreq: %d', options.bandPass)];
                end
                title(titleText)
            end
        end


        function plot_trials(obj, animalId, recName, channel, pad)
            arguments
                obj predictionBreak
                animalId char
                recName char
                channel double = nan
                pad double = 1000  % extra time before and after in msec
            end
            [arena, T, defaultCh] = obj.load_rec(animalId, recName);
            if isnan(channel), channel = defaultCh; end
            [tEnter, ~] = obj.get_event_in_OE_time(arena, T.enter_time, 'enter_time', animalId, recName);
            [tEscape, ~] = obj.get_event_in_OE_time(arena, T.escape_time, 'escape_time', animalId, recName);
            [tExit, ~] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time', animalId, recName);

            cols = 3;
            rows =  ceil(numel(tEnter)/cols);
            figure();
            tlo = tiledlayout(rows, cols, 'TileSpacing','compact', 'Padding','compact');
            title(tlo, sprintf('%s, %s, channel=%d', animalId, recName, channel));
            for i=1:numel(tEnter)
                nexttile;
                duration = tExit(i) - tEnter(i);
                [V, t] = obj.wa.currentDataObj.getData(channel, tEnter(i)-pad, duration+2*pad);
                V = obj.wa.filt.FL.getFilteredData(V);
                plot(t, squeeze(V));
                xline(pad, '--r')
                xline(pad+tEscape(i)-tEnter(i), '--g')
                xline(pad+tExit(i)-tEnter(i), '--r')
            end
        end


        function display_events(obj, animalId, recName, cam_name, eventName)
            [arena, T, ~] = load_rec(obj, animalId, recName);
            T = T(~isnat(T.(eventName)), :);
            [~, event_frames] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, recName);
            event_labels = cell(1, length(event_frames));
            tagCol = sprintf('%s_tags', eventName);
            for i=1:height(T)
                lbl = sprintf('%d trial%d', event_frames(i), T.trial_id(i));
                if strcmp(T.(tagCol){i}, 'not engaged')
                    lbl = [lbl ' (not_eng)'];
                end
                if strcmp(T.(tagCol){i}, 'strike')
                    lbl = [lbl ' (strk)'];
                end
                event_labels{i} = lbl;
            end
            if isempty(event_frames)
                fprintf('No events were found for %s\n', recName)
                return
            end
            video_path = find_video_files(obj.wa.currentExpFolder, cam_name);
            if isnan(video_path)
                return; % Exit if no video file found
            end
            title = sprintf('%s %s', animalId, recName);
            frame_navigator(video_path, event_frames(1), event_frames, event_labels, title)
        end
        

        function play_trial(obj, animalId, recName, trialId, cam_name)
            [arena, T, ~] = load_rec(obj, animalId, recName);
            T = T(T.trial_id==trialId, :);
            [~, enter_frame] = obj.get_event_in_OE_time(arena, T.enter_time, 'enter_time', '', recName);
            [~, exit_frame] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time', '', recName);
            video_path = find_video_files(obj.wa.currentExpFolder, cam_name);
            if isnan(video_path)
                return; % Exit if no video file found
            end
            play_video_segment(video_path,enter_frame,exit_frame-enter_frame)
        end


        function [event_oe_time, event_frame] = get_event_in_OE_time(obj, arena, eventTime, eventName, animalId, recName)
            arguments
                obj predictionBreak
                arena struct
                eventTime datetime
                eventName char
                animalId char = ''
                recName char = ''
            end
            dt = datetime(eventTime, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'Asia/Jerusalem');
            ts = posixtime(dt);
            event_frame = obj.wa.getVideoFrames(arena.videoFrames, ts);
            event_frame = event_frame - arena.IRFrames(1);
            
            if ~isempty(animalId) & ~isempty(recName)
                Gp = obj.frameGap;
                Gp = Gp(strcmp(Gp.animal_id, animalId) & strcmp(Gp.rec_name, recName), :);
                if ~isempty(Gp) & ~isnan(Gp.(eventName))
                    event_frame = event_frame - Gp.(eventName);
                end
            end
            event_oe_time = arena.oeCamTrigs(event_frame);
        end


        function [st, clu, cluster_annot, arena, defaultCh] = get_spikes_data(obj, animalId, recName, isRemoveNoise)
            arguments
                obj predictionBreak
                animalId char
                recName char
                isRemoveNoise logical = true
            end
            [arena, ~, defaultCh] = obj.load_rec(animalId, recName);
            results_dir = dir(fullfile(obj.wa.currentExpFolder, 'spikeSorting', 'kilosort4'));
            if isempty(results_dir)
                fprintf('No spike sorting results folder found in %s\n', obj.wa.currentExpFolder)
                return
            end
            results_dir = results_dir(1).folder;
            st  = readNPY(fullfile(results_dir,'spike_times.npy'));       % samples
            clu = readNPY(fullfile(results_dir,'spike_clusters.npy'));
            cluster_annot = readtable(fullfile(results_dir,'cluster_group.tsv'),'FileType','text');

            % remove all noise clusters
            if isRemoveNoise
                good_clu = cluster_annot.cluster_id(~strcmp(cluster_annot.group, 'noise'));
                st = st(ismember(clu, good_clu));
                clu = clu(ismember(clu, good_clu));
            end
        end
        

        function plot_spikes_per_event(obj, animalId, recName, eventName, options)
            arguments
                obj predictionBreak
                animalId char
                recName char
                eventName char
                options.is_engaged logical = false
                options.sec_before double = 3
                options.sec_after double = 3
                options.binw double = 100
                options.ymax double = 100
                options.is_remove_noise = true
            end
            tiledlayout('flow','TileSpacing','compact','Padding','compact');
            [st, clu, cluster_annot, arena, defaultCh] = obj.get_spikes_data(animalId, recName, options.is_remove_noise);
            
            fs = obj.wa.currentDataObj.samplingFrequency(defaultCh);
            t_ms = 1000 * double(st) / fs;

            T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, recName);
            [tEvent, ~] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, recName);
            if strcmp(eventName, 'escape_time')
                [tEventExit, ~] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time', animalId, recName);
            end
            
            tStart = tEvent - (options.sec_before * 1000);
            tEnd = tEvent + (options.sec_after * 1000);

            nI = numel(tStart);
            counts = cell(nI,1);
            edges  = cell(nI,1);
            for k = 1:nI
                a = tStart(k); b = tEnd(k);
                if a > b, [a,b] = deal(b,a); end
                % Bin edges covering [a,b] with width binw (last bin right-inclusive)
                nBins = max(1, ceil((b - a)/options.binw));
                e = a + (0:nBins)*options.binw;
                if e(end) < b, e(end+1) = b; end
                % Data inside interval
                x = t_ms(t_ms >= a & t_ms <= b);
                % Counts
                counts{k} = histcounts(x, e);
                edges{k}  = e;
                % Plot
                nexttile;
                histogram((x-tEvent(k))/1000, 'BinEdges', (e-tEvent(k))/1000);  % , 'Normalization', 'probability'
                xlabel(sprintf('Time (sec)')); ylabel('Count');
                title(sprintf('Interval %d: [%g, %g], bin=%.3gms', k, a, b, options.binw));
                ylim([0 options.ymax])
                xline(0, '--k')
                if strcmp(eventName, 'escape_time')
                    xline((tEventExit(k)-tEvent(k))/1000, '--r')
                end
                grid on
            end
        end


        function plot_avg_event_spikes_per_cluster(obj, animalId, recName, eventName, options)
            arguments
                obj predictionBreak
                animalId char
                recName char
                eventName char
                options.is_engaged logical = false
                options.sec_before double = 3
                options.sec_after double = 3
                options.binw double = 100 % in milliseconds
                options.is_remove_noise = true
            end
            [st, clu, cluster_annot, arena, defaultCh] = obj.get_spikes_data(animalId, recName, options.is_remove_noise);
            uniq_clusters = sort(unique(clu));
            
            fs = obj.wa.currentDataObj.samplingFrequency(defaultCh);
            t_ms = 1000 * double(st) / fs;
            T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, recName);
            [tEvent, ~] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, recName);
            if strcmp(eventName, 'escape_time')
                [tEventExit, ~] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time', animalId, recName);
            end
            tStart = tEvent - (options.sec_before * 1000);
            tEnd = tEvent + (options.sec_after * 1000);
            nI = numel(tStart);
            
            figure();
            tdl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
            title(tdl, sprintf('%s, %s, %s', animalId, recName, eventName));
            for i = 1:length(uniq_clusters)
                clusterId = uniq_clusters(i);
                counts = cell(nI,1);
                edges  = cell(nI,1);
                for k = 1:nI
                    a = tStart(k); b = tEnd(k);
                    if a > b, [a,b] = deal(b,a); end
                    % Bin edges covering [a,b] with width binw (last bin right-inclusive)
                    nBins = max(1, ceil((b - a)/options.binw));
                    e = a + (0:nBins)*options.binw;
                    if e(end) < b, e(end+1) = b; end
                    % Data inside interval
                    t_ms_clu = t_ms(clu==clusterId);
                    x = t_ms_clu(t_ms_clu >= a & t_ms_clu <= b);
                    % Counts
                    counts{k} = histcounts(x, e);
                    edges{k}  = (e - tEvent(k))/1000;
                end

                % Plot
                nexttile;
                counts = vertcat(counts{:});
                rate = counts / (options.binw*1000);
                % t = (0:nBins-1)*options.binW + options.binW/2;
                t = edges{1}(1:end-1);
                hold on
                plot(t, rate.', 'Color', [0.8 0.8 0.8]);          % each segment (light)
                m  = mean(rate, 1, 'omitnan');                    % 1 × nBins
                s  = std(rate, 0, 1, 'omitnan') / sqrt(size(rate,1));  % SEM
                % Shaded SEM
                fill([t fliplr(t)], [m-s fliplr(m+s)], [0.3 0.6 0.9], ...
                     'FaceAlpha', 0.25, 'EdgeColor', 'none');
                % Mean curve
                plot(t, m, 'Color', [0.1 0.35 0.7], 'LineWidth', 2);
                xlabel('Time within segment (s)');
                ylabel('Spike rate (Hz)');
                box off
                title(sprintf('Cluster %d (%s)', clusterId, cluster_annot.group{cluster_annot.cluster_id==clusterId}));
                xline(0, '--k')
            end
        end
        
        
        function [recNames] = get_recs_for_animal(obj, animal_id, eventName, movementType)
            recNames = unique(obj.mainT.rec_name(strcmp(obj.mainT.animal_id, animal_id)&...
                ~isnat(obj.mainT.(eventName)) & strcmp(obj.mainT.movement_type, movementType)));
        end
        

        function T = get_table_for_event_name(obj, animalId, eventName, is_engaged, rec_name, movementType, no_strikes)
            arguments
                obj predictionBreak
                animalId char
                eventName char
                is_engaged logical
                rec_name char = ''
                movementType char = ''
                no_strikes logical = true
            end
            T = obj.mainT(strcmp(obj.mainT.animal_id, animalId), :);
            if ~isempty(movementType)
                T = T(strcmp(T.movement_type, movementType), :);
            end 
            if ~isempty(rec_name)
                T = T(strcmp(T.rec_name, rec_name), :);
            end
            T = T(~isnat(T.(eventName)), :);
            if strcmp(eventName, 'exit_time')
                T = T(~isnat(T.escape_time), :);
            end
            tagCol = sprintf('%s_tags', eventName);
            if is_engaged
                % T = T(strcmp(T.(sprintf('is_%s_engaged', erase(eventName, '_time'))), 'True'), :);
                T = T(~strcmp(T.(tagCol), 'not engaged'), :);
            end
            if no_strikes
                % T = T(T.n_escape_strikes==0 | isnan(T.n_escape_strikes), :);
                T = T(~strcmp(T.(tagCol), 'strike'), :);
            end
        end

        
        function rec_names = get_relevant_rec_names(obj, animalId, eventName, movementType, options)
            arguments
                obj predictionBreak
                animalId char
                eventName char
                movementType char
                options.is_engaged logical = true
                options.no_strikes logical = true
                options.min_trials double = nan
            end
            T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, '', movementType, options.no_strikes);
            if ~isnan(options.min_trials)
                summaryT = groupsummary(T, "rec_name");
                rec_names = summaryT.rec_name(summaryT.GroupCount >= options.min_trials);
            else
                rec_names = unique(T.rec_name);
            end
            % sort rec names by hunter number
            nums = cellfun(@(s) sscanf(s,'Hunter%d'), rec_names);
            [~, idx] = sort(nums);
            rec_names = rec_names(idx);
        end


        function plot_all_events_with_images(obj, animalId, eventName, movementType, options)
            arguments
                obj predictionBreak
                animalId char
                eventName char {mustBeMember(eventName, ['enter_time','escape_time', 'exit_time', 'flip_time'])}
                movementType char
                options.channel double = nan
                options.time_before double = 2 % time in sec
                options.time_after double = 2 % time in sec
                options.is_engaged logical = true
                options.no_strikes logical = true
                options.max_cols_per_page double = 5;
                options.cam_name char = 'right'
                options.gammaLow double = 40
                options.gammaHigh double = 100
            end
            isEscape = strcmp(eventName, 'escape_time');
            T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, '', ...
                                             movementType, options.no_strikes);
            recNames = unique(T.rec_name);
            nTrials = height(T);
            nrows = 3;
            totalCols = ceil(nTrials/nrows);
            
            fig = figure('PaperUnits','centimeters','PaperPosition',[0 0 21 29]);
            tlo = tiledlayout(nrows, options.max_cols_per_page, 'TileSpacing', 'compact', 'Padding', 'normal','TileIndexing','columnmajor');
            currentCol = 1;
            currentPage = 1;
            k = 1;
            for i=1:numel(recNames)
                try
                    [arena, ~, channel] = obj.load_rec(animalId, recNames{i});
                    Trec = T(strcmp(T.rec_name, recNames{i}), :);
                    [tEvent, vidFrames] = obj.get_event_in_OE_time(arena, Trec.(eventName), eventName, animalId, recNames{i});
                    if isEscape
                        [tExit, ~] = obj.get_event_in_OE_time(arena, Trec.exit_time, 'exit_time', animalId, recNames{i});
                    end
                    [V, ~] = obj.wa.currentDataObj.getData(channel, tEvent-(options.time_before*1000), (options.time_before+options.time_after)*1000);
                    V = squeeze(V);
                    if size(V, 2) == 1 % squeeze flip order if only one segment
                        V = V';
                    end
                    % load video
                    video_path = find_video_files(obj.wa.currentExpFolder, options.cam_name);
                    vid = VideoReader(video_path);

                    for j=1:length(tEvent)
                        % plot spectrogram
                        nexttile;
                        try
                            fs = obj.wa.currentDataObj.samplingFrequency(channel);
                            S = obj.run_spectrogram(V(j, :), fs, time_before=options.time_before, ...
                                time_after=options.time_after, clim=[-30,20],...
                                ax=true, freqRange=[1 150]);
                            titleText = sprintf('%s, trial%d',recNames{i},Trec.trial_id(j));
                            tagCol = sprintf('%s_tags', eventName);
                            if ~options.is_engaged && ~strcmp(Trec.(tagCol){j}, 'not engaged')
                                titleText = [titleText ' (engaged)'];
                            end
                            if ~options.no_strikes && ~strcmp(Trec.(tagCol){j}, 'strike')
                                titleText = [titleText ' (strike)'];
                            end
                            title(titleText)
                            if isEscape
                                xline((tExit(j)-tEvent(j))/1000, '--r')
                            end
                        catch ME
                            fprintf('Error in plotting spectrogram rec %s, trial %d: %s\n', recNames{i}, Trec.trial_id(j), ME.message)
                        end
                        % plot image
                        nexttile;
                        try
                            vid.CurrentTime = (vidFrames(j)+3) / vid.FrameRate;
                            img = readFrame(vid);
                            imshow(img)
                            axis off
                        catch ME
                            fprintf('Error in plotting image rec %s, trial %d: %s\n', recNames{i}, Trec.trial_id(j), ME.message)
                        end

                        % plot beta2gamma
                        nexttile;
                        try
                            gamma = mean(S.P_event((options.gammaLow<=S.f)&(S.f<=options.gammaHigh),:), 1);
                            plot(S.t, gamma)
                            xline(0, '--k')
                            if isEscape
                                t_exit_sec = (tExit(j)-tEvent(j))/1000;
                                if t_exit_sec <= options.time_after
                                    xline(t_exit_sec, '--r');
                                end
                            end
                            set(gca,'Color','w')
                        catch ME
                            fprintf('Error in plotting gamma; rec %s, trial %d: %s\n', recNames{i}, Trec.trial_id(j), ME.message)
                        end

                        currentCol = currentCol + 1;
                        if (currentCol > options.max_cols_per_page) || ((i==numel(recNames)) && (j==length(tEvent)))
                            exportgraphics(fig, sprintf('%s/%s_%s_%s.pdf', obj.outputDir, animalId, eventName, movementType), 'ContentType','vector', 'Append', currentPage>1)
                            close(fig)
                            fig = figure('PaperUnits','centimeters','PaperPosition',[0 0 25 29]);
                            tlo = tiledlayout(nrows, options.max_cols_per_page, 'TileSpacing', 'compact', 'Padding', 'normal', 'TileIndexing','columnmajor');
                            currentPage = currentPage + 1;
                            currentCol = 1;
                        end
                    end
                catch ME
                    fprintf('Error rec %s: %s\n', recNames{i}, ME.message)
                    % rethrow(ME)
                end
            end
        end

        
        function print_animals(obj)
            animals = unique(obj.mainT.animal_id);
            for i=1:numel(animals)
                fprintf('%s:\n', animals{i})
                T = obj.mainT(strcmp(obj.mainT.animal_id, animals{i}), :);
                movements = unique(T.movement_type);
                for j=1:numel(movements)
                    T1 = T(strcmp(T.movement_type, movements{j}), :);
                    fprintf('  %s:\n', movements{j})
                    escapeT = T1(~isnat(T1.escape_time), :);
                    tagCol = escapeT.escape_time_tags;
                    fprintf('    escape: %d/%d\n', sum(~strcmp(tagCol, 'not engaged') & ~strcmp(tagCol, 'strike')), height(escapeT));
                    
                    flipT = T1(~isnat(T1.flip_time), :);
                    tagCol = flipT.flip_time_tags;
                    if ~isempty(flipT)
                        fprintf('    flip: %d/%d\n', sum(~strcmp(tagCol,'not engaged') & ~strcmp(tagCol, 'strike')), height(flipT));
                    end
                end
            end
        end
        
        
        function cacheFile = get_cache_file_name(obj, name, options, options2drop, varargin)
            options2drop = intersect(fieldnames(options), cellstr(string(options2drop)));
            options = rmfield(options, options2drop);
            vals = struct2cell(options);
            % Convert each value to string
            valsStr = cellfun(@num2str, vals, 'UniformOutput', false);
            outStr = strjoin(valsStr, '_');
            cacheFile = sprintf('%s/%s_%s_%s.mat', obj.cacheDir, name, strjoin(varargin, '_'), outStr);
        end


        function plot_video_frames(obj, cam_name, frame_ids, n_frames)
            video_path = find_video_files(obj.wa.currentExpFolder, cam_name);
            if isnan(video_path)
                return; % Exit if no video file found
            end
            v = VideoReader(video_path);
            num_cols = min(numel(frame_ids), 7);
            % Create figure
            figure('Units', 'normalized', 'Position', [0, 0, 1, 1])
            t = tiledlayout(n_frames, num_cols, 'TileSpacing', 'none', 'Padding', 'none');
            for i = 1:n_frames
                for j = 1:num_cols
                    frame_num = frame_ids(j) + (i-1); % get frame index
                    v.CurrentTime = (frame_num-1) / v.FrameRate;
                    img = readFrame(v);
        
                    nexttile
                    imshow(img)
                    axis off
                end
            end
        end
        

        function [t_active,lags] = get_screen_diode_activations(obj, arena, T, isPlot)
            arguments
                obj predictionBreak
                arena struct
                T table
                isPlot logical = true
            end
            t_active = nan(height(T), 2);
            lags = nan(height(T), 2);
            if isempty(T)
                fprintf('Trials table is empty')
                return
            end
            [tStarts, ~] = obj.get_event_in_OE_time(arena, T.enter_time, 'enter_time');
            [tEnds, ~] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time');
            oe = obj.wa.currentDataObj;
            dt = 3000;
            % trial duration shouldn't be more than 43 seconds
            trialDuration = min(max(tEnds-tStarts), 43e3);
            [V_diode, t_diode] = oe.getAnalogData(1, tStarts-dt, 2*dt+trialDuration);
            V_diode = squeeze(V_diode);
            if isscalar(tStarts) % squeeze flip order if only one segment
                V_diode = V_diode';
            end
            if isPlot
                cols = 2;
                rows = ceil(height(V_diode)/cols);
                figure;
                tiledlayout(rows, cols, "TileSpacing","tight","Padding","tight");
            end
            for i=1:height(V_diode)
                v = V_diode(i, :);
                th = mean(v) - 6*std(v);
                below = v < th;
                d = diff([0; below(:); 0]);  % pad with 0s at edges
                starts = find(d == 1);
                ends = find(d == -1) - 1;
                lengths = ends - starts + 1;
                validStarts = starts(lengths >= 100);
        
                if isPlot
                    nexttile();
                    plot(t_diode, v)
                    xline(dt, '--g')
                    xline(dt + (tEnds(i)-tStarts(i)), '--g')
                    yline(th, '--r')
                    hold on
                    scatter(t_diode(validStarts), v(validStarts), 100, 'filled', 'MarkerFaceColor', 'm')
                    % print title
                    cond = length(validStarts) > 1;
                    duration = (t_diode(validStarts(end)) - t_diode(validStarts(1)))*cond + 0*~cond;
                    title(sprintf('%.1fsec', duration/1000));
                end
        
                % calcualte the OE lag
                lags(i, 1) = dt - t_diode(validStarts(1));
                lags(i, 2) = dt + (tEnds(i)-tStarts(i)) - t_diode(validStarts(end));
        
                t = t_diode + (tStarts(i)-dt);
                t_active(i, 1) = t(validStarts(1));
                t_active(i, 2) = t(validStarts(end));
            end
        end
    end
end



        % 
        % 
        % function [common_t, all_b2g] = get_beta2Gamma(obj, animalId, recName, eventName, options)
        %     arguments
        %         obj predictionBreak
        %         animalId char
        %         recName char
        %         eventName char
        %         options.is_engaged logical = false
        %         options.time_before double = 10
        %         options.time_after double = 5
        %         options.gammaLow double = 40
        %         options.gammaHigh double = 100
        %         options.betaLow double = 8
        %         options.betaHigh double = 20
        %         options.movWin double = 1000
        %         options.movOLWin double = 800
        %         options.overwrite double = 0
        %         options.is_plot logical = true
        %         options.use_spectrogram logical = true
        %     end
        %     [arena, ~, channel] = obj.load_rec(animalId, recName);            
        %     T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, recName);
        %     [tEvent, ~] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, recName);
        %     isEscape = strcmp(eventName, 'escape_time');
        %     if isEscape
        %         [tExit, ~] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time', animalId, recName);
        %     end
        % 
        %     step = (options.movWin - options.movOLWin) / 1000;
        %     common_t = -options.time_before:step:options.time_after;
        % 
        %     if options.use_spectrogram
        %         all_b2g = obj.get_beta2gamma_by_spectrogram(tEvent, channel, common_t, time_before=options.time_before, ...
        %             time_after=options.time_after, gammaLow=options.gammaLow, gammaHigh=options.gammaHigh, ...
        %             betaLow=options.betaLow, betaHigh=options.betaHigh, movOLWin=options.movOLWin, movWin=options.movWin);
        %     else
        %         all_b2g = obj.get_beta2gamma_by_wakeAnalysis(tEvent, common_t, time_before=options.time_before, ...
        %             time_after=options.time_after, gammaLow=options.gammaLow, gammaHigh=options.gammaHigh, ...
        %             betaLow=options.betaLow, betaHigh=options.betaHigh, movOLWin=options.movOLWin, movWin=options.movWin, ...
        %             overwrite=options.overwrite);
        %     end
        % 
        %     if options.is_plot
        %         figure();
        %         tdl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
        %         title(tdl, sprintf('%s - %s - %s\nbeta:%d-%d, gamma:%d-%d', animalId, recName, eventName, options.betaLow,options.betaHigh, ...
        %             options.gammaLow,options.gammaHigh),'Interpreter','none')
        %         for i=1:length(tEvent)
        %             nexttile;
        %             plot(common_t, all_b2g(i, :))
        %             xline(0, '--k')
        %             if isEscape
        %                 xline((tExit(i)-tEvent(i))/1000, '--r')
        %             end
        %             title(sprintf('trial%d', T.trial_id(i)))
        %         end
        %         ax = findobj(tdl, 'Type','axes');
        %         linkaxes(ax,'y');
        %     end
        % end
        % 
        % 
        % function all_b2g = get_beta2gamma_by_spectrogram(obj, tEvent, channel, common_t, options)
        %     arguments
        %         obj predictionBreak
        %         tEvent double
        %         channel double
        %         common_t double
        %         options.time_before double = 10
        %         options.time_after double = 5
        %         options.gammaLow double = 40
        %         options.gammaHigh double = 100
        %         options.betaLow double = 8
        %         options.betaHigh double = 20
        %         options.movWin double = 1000
        %         options.movOLWin double = 800
        %     end
        %     fs = obj.wa.currentDataObj.samplingFrequency(channel);
        %     duration = options.time_before + options.time_after;
        %     all_b2g = nan(length(tEvent), length(common_t));
        %     for i=1:length(tEvent)
        %         [V, ~] = obj.wa.currentDataObj.getData(channel, tEvent(i)-(options.time_before*1000), duration*1000);
        %         V = squeeze(V);
        %         if size(V, 2) == 1 % squeeze flip order if only one segment
        %             V = V';
        %         end
        %         S = event_aligned_spectrogram(V, fs, ...
        %             'FreqRange', [1 100], ...
        %             'WindowSec', options.movWin/1000, ...
        %             'StepSec', (options.movWin-options.movOLWin)/1000, ...
        %             'PreWindow', [-options.time_before 0], ...
        %             'FullWindow',[-options.time_before options.time_after], ...
        %             'Plot', false);
        %         beta = mean(S.P_event((options.betaLow<=S.f)&(S.f<=options.betaHigh),:), 1);
        %         gamma = mean(S.P_event((options.gammaLow<=S.f)&(S.f<=options.gammaHigh),:), 1);
        %         all_b2g(i, :) = interp1(S.t,beta./gamma,common_t,'spline');
        %     end
        % end
        % 
        % 
        % function all_b2g = get_beta2gamma_by_wakeAnalysis(obj, tEvent, common_t, options)
        %     arguments
        %         obj predictionBreak
        %         tEvent double
        %         common_t double
        %         options.time_before double = 10
        %         options.time_after double = 5
        %         options.gammaLow double = 40
        %         options.gammaHigh double = 100
        %         options.betaLow double = 8
        %         options.betaHigh double = 20
        %         options.movWin double = 1000
        %         options.movOLWin double = 800
        %         options.overwrite double = 0
        %     end
        %     b2g = obj.wa.getBeta2GammaRatio('gammaLow', options.gammaLow, 'gammaHigh', options.gammaHigh, 'betaLow', options.betaLow, ...
        %                                     'betaHigh', options.betaHigh, 'movWin', options.movWin,'movOLWin', options.movOLWin, ...
        %                                     'overwrite', options.overwrite);
        %     all_b2g = nan(length(tEvent), length(common_t));
        %     for i=1:length(tEvent)
        %         idx = (tEvent(i) - (options.time_before*1000) <= b2g.t_ms) & (b2g.t_ms <= tEvent(i) + (options.time_after*1000));
        %         t = (b2g.t_ms(idx)-tEvent(i))/1000;
        %         y = b2g.beta2gammaRatio(idx);
        %         all_b2g(i, :) = interp1(t,y,common_t,'spline');
        %     end
        % end
        % 
        % 
        % function plot_all_beta_to_gamma(obj, animalId, eventName, movementType, options)
        %     arguments
        %         obj predictionBreak
        %         animalId char
        %         eventName char
        %         movementType char = 'circle'
        %         options.is_engaged logical = true
        %         options.time_before double = 10
        %         options.time_after double = 5
        %         options.no_strikes logical = false
        %         options.all_together logical = false
        %         options.min_trials double = nan
        %         options.use_spectrogram logical = true
        %         options.ax logical = false
        %         options.gammaLow double = 40
        %         options.gammaHigh double = 100
        %         options.betaLow double = 8
        %         options.betaHigh double = 20
        %         options.movWin double = 1000
        %         options.movOLWin double = 800
        %         options.overwrite double = 0
        %     end
        %     fprintf('>> start plotting all beta2gamma for %s, %s, %s\n', animalId, movementType, eventName);
        %     rec_names = obj.get_relevant_rec_names(animalId, eventName, movementType, no_strikes=options.no_strikes, ...
        %                                            is_engaged=options.is_engaged, min_trials=options.min_trials);
        %     if ~options.ax
        %         figure();
        %         tdl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
        %         title(tdl, sprintf('%s - %s', animalId, eventName),'Interpreter','none')
        %     end
        % 
        %     cacheFile = obj.get_cache_file_name('b2g', options, animalId, eventName, movementType);
        %     if isfile(cacheFile) && options.all_together
        %         load(cacheFile, 'res', 'common_t'); 
        %     else
        %         res = [];
        %         for i=1:numel(rec_names)
        %             if ~mod(i, 25)
        %                 fclose('all');
        %             end
        %             try
        %                 [common_t, all_b2g] = obj.get_beta2Gamma(animalId, rec_names{i}, eventName, is_plot=false, is_engaged=options.is_engaged, ...
        %                                                          use_spectrogram=options.use_spectrogram, gammaLow=options.gammaLow, ...
        %                                                          gammaHigh=options.gammaHigh, betaLow=options.betaLow, betaHigh=options.betaHigh, ...
        %                                                          movOLWin=options.movOLWin, movWin=options.movWin, ...
        %                                                          overwrite=options.overwrite);
        %                 res = [res; all_b2g];
        %                 if ~options.all_together
        %                     nexttile;
        %                     plot(common_t, mean(all_b2g, 1))
        %                     title(sprintf('%s', rec_names{i}))
        %                 end
        %             catch ME
        %                 fprintf('Error in rec %s (%d/%d): %s\n', rec_names{i}, i, numel(rec_names), ME.message)
        %             end
        %         end
        %         % remove rows with all NaN
        %         res(all(isnan(res),2),:) = [];
        %         save(cacheFile, 'res', 'common_t'); 
        %     end
        %     if options.all_together
        %         plot(common_t, mean(res, 1))
        %         xline(0, 'k')
        %         title_text = sprintf('%s - %s - %s - n_trials=%d', animalId, movementType, eventName, size(res,1));
        %         if options.is_engaged
        %             title_text = [title_text ' (engaged)'];
        %         end
        %         title_text = [title_text sprintf('\nbeta:%d-%d, gamma:%d-%d', options.betaLow,options.betaHigh, options.gammaLow,options.gammaHigh)];
        %         if options.ax
        %             title(title_text, 'Interpreter','none');
        %         else
        %             title(tdl, title_text, 'Interpreter','none');
        %         end
        %     end
        % end
        % 
        % 
        % function plot_b2g_and_spect_for_entire_trial(obj, animalId, recName, trialId, eventName, options)
        %     arguments
        %         obj predictionBreak
        %         animalId char
        %         recName char
        %         trialId double
        %         eventName char
        %         options.is_engaged logical = false
        %         options.extra_time double = 2
        %         options.gammaLow double = 60
        %         options.gammaHigh double = 100
        %         options.betaLow double = 8
        %         options.betaHigh double = 20
        %         options.movWin double = 1000
        %         options.movOLWin double = 800
        %         options.overwrite double = 0
        %         options.is_plot logical = true
        %         options.time_around_event double = nan
        %     end
        %     [arena, ~, channel] = obj.load_rec(animalId, recName);
        %     % b2g = obj.wa.getBeta2GammaRatio('gammaLow', options.gammaLow, 'gammaHigh', options.gammaHigh, 'betaLow', options.betaLow, ...
        %     %                                 'betaHigh', options.betaHigh, 'movWin', options.movWin,'movOLWin', options.movOLWin, ...
        %     %                                 'overwrite', options.overwrite);
        %     T = obj.get_table_for_event_name(animalId, eventName, options.is_engaged, recName);
        %     T = T(T.trial_id==trialId,:);
        %     [tEvent, ~] = obj.get_event_in_OE_time(arena, T.(eventName), eventName, animalId, recName);
        %     [tExit, ~] = obj.get_event_in_OE_time(arena, T.exit_time, 'exit_time', animalId, recName);
        % 
        %     if isnan(options.time_around_event)
        %         [tStart, ~] = obj.get_event_in_OE_time(arena, T.enter_time, 'enter_time', animalId, recName);
        %         duration = (2*options.extra_time+(tExit-tStart)/1000);
        %         startTime = tStart-(options.extra_time*1000);
        %         time_before = options.extra_time;
        %         time_after = duration-options.extra_time;
        %     else
        %         duration = 2*options.time_around_event;
        %         startTime = tEvent - (options.time_around_event*1000);
        %         tStart = tEvent;
        %         time_before = options.time_around_event;
        %         time_after = options.time_around_event;
        %     end
        % 
        %     [V, t_v] = obj.wa.currentDataObj.getData(channel, startTime, duration*1000);
        %     V = squeeze(V);
        %     if size(V, 2) == 1 % squeeze flip order if only one segment
        %         V = V';
        %     end
        % 
        %     figure();
        %     tld = tiledlayout(1,3);
        %     title(tld, sprintf('%s - %s - trial%d', animalId, recName, trialId), 'Interpreter','none')
        %     % spectrogtam
        %     nexttile;
        %     S = obj.run_spectrogram(V, channel, animalId, recName, eventName, WindowSec=options.movWin/1000, ...
        %                                StepSec=(options.movWin-options.movOLWin)/1000, ax=true, time_before=time_before, ...
        %                                time_after=time_after, clim=[-30 20]);
        %     xline((tExit-tStart)/1000, 'k')
        %     xline((tEvent-tStart)/1000, '--r')
        %     title('spectrogram')
        % 
        %     % beta2gamma from spectrogram
        %     nexttile;
        %     beta = mean(S.P_event((options.betaLow<=S.f)&(S.f<=options.betaHigh),:), 1);
        %     gamma = mean(S.P_event((options.gammaLow<=S.f)&(S.f<=options.gammaHigh),:), 1);
        %     plot(S.t, beta./gamma)
        %     xline(0, 'k')
        %     xline((tExit-tStart)/1000, 'k')
        %     xline((tEvent-tStart)/1000, '--r')
        % 
        %     % beta2gamma from code
        %     % nexttile;
        %     % idx = (tStart - (options.extra_time*1000) <= b2g.t_ms) & (b2g.t_ms <= tExit + (options.extra_time*1000));
        %     % t = (b2g.t_ms(idx)-tStart)/1000;
        %     % y = b2g.beta2gammaRatio(idx);
        %     % plot(t, y)
        %     % xline(0, 'k')
        %     % xline((tExit-tStart)/1000, 'k')
        %     % xline((tEvent-tStart)/1000, '--r')
        %     % title(sprintf('beta:%d-%d, gamma:%d-%d', options.betaLow,options.betaHigh, ...
        %     %         options.gammaLow,options.gammaHigh),'Interpreter','none')
        % end

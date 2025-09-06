%% print animal data

pb = predictionBreak;
pb.print_animals()

%% dispaly events from videos

animal_id = 'PV153';
event_name = 'escape_time';
pb = predictionBreak;
recNames = pb.get_recs_for_animal(animal_id, event_name, 'low_horizontal_noise');
pb.display_events(animal_id, recNames{1}, 'right', event_name)

% PV106, Hunter40 study first escape - strike

%% events spectrograms

animal_id = 'PV153';
pb = predictionBreak;
% recNames = pb.get_recs_for_animal(animal_id);
% pb.calc_event_spectrogram_for_rec(animal_id,'Hunter48','escape_time',time_before=2,time_after=2,is_engaged=false,freqRange=[1 300]);

pb.plot_all_event_spectrogram(animal_id, 'circle', 'escape_time', min_trials=1, time_before=1, time_after=1, ...
                              all_together=true, is_engaged=true, freqRange=[1 150], is_cache=true, isPlot=true);

%% play trial
pb = predictionBreak;
pb.play_trial('PV126', 'Hunter21', 15, 'right')

%% plot segments

animal_id = 'PV106';
pb = predictionBreak;
% pb.plot_event_segments(animal_id, 'Hunter17', 'escape_time', 12, bandRange=[30, 100]) % best circle escape
% pb.plot_event_segments(animal_id, 'Hunter17', 'escape_time', 6)  % , highPass=1, lowPass=1000
% pb.plot_event_segments(animal_id, 'Hunter7', 'escape_time', 3)

% pb.plot_event_segments(animal_id, 'Hunter17', 'flip_time', 15)

pb.plot_event_segments('PV157', 'Hunter72', 'escape_time', 14, bandRange=[40, 100])
% pb.plot_event_segments('PV157', 'Hunter73', 'escape_time', 16, bandRange=[30, 150])

%% plot all trials with images

pb = predictionBreak;
% pb.plot_all_events_with_images('PV153', 'escape_time', 'circle', cam_name='right')
% pb.plot_all_events_with_images('PV153', 'escape_time', 'low_horizontal_noise', cam_name='right')
pb.plot_all_events_with_images('PV153', 'flip_time', 'circle', cam_name='right')

% close all
% 
% animals = {'PV106', 'PV157', 'PV126'};
% for i=1:numel(animals)
%     pb = predictionBreak;
%     pb.plot_all_events_with_images(animals{i}, 'escape_time', 'circle', cam_name='right')
%     close all
% end
% 
% animals = {'PV157', 'PV126'};
% for i=1:numel(animals)
%     pb = predictionBreak;
%     pb.plot_all_events_with_images(animals{i}, 'escape_time', 'low_horizontal_noise', cam_name='right')
%     close all
% end
% 
% animals = {'PV106', 'PV143'};
% for i=1:numel(animals)
%     pb = predictionBreak;
%     pb.plot_all_events_with_images(animals{i}, 'flip_time', 'circle', cam_name='right')
%     close all
% end


%% plot all spectrograms on one figure

figure();
tdl = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
animals = {'PV157', 'PV106', 'PV126'};
for i=1:numel(animals)
    nexttile;
    pb = predictionBreak;
    pb.plot_all_event_spectrogram(animals{i}, 'circle', 'escape_time', min_trials=1, time_before=1, time_after=1, all_together=true, is_engaged=true, freqRange=[1 150], ax=true);
end
animals = {'PV157', 'PV162'};
for i=1:numel(animals)
    nexttile;
    pb = predictionBreak;
    pb.plot_all_event_spectrogram(animals{i}, 'low_horizontal_noise', 'escape_time', min_trials=1, time_before=1, time_after=1, all_together=true, is_engaged=true, freqRange=[1 150], ax=true);
end
nexttile;
animals = {'PV106', 'PV143'};
for i=1:numel(animals)
    nexttile;
    pb = predictionBreak;
    pb.plot_all_event_spectrogram(animals{i}, 'circle', 'flip_time', min_trials=1, time_before=1, time_after=1, all_together=true, is_engaged=true, freqRange=[1 150], ax=true);
end

%% spike sorting example

pb = predictionBreak;
% pb.plot_avg_event_spikes_per_cluster('PV126', 'Hunter23', 'escape_time', is_remove_noise=false)

% pb.plot_spikes_per_event('PV157', 'Hunter42', 'escape_time')
% pb.plot_avg_event_spikes_per_cluster('PV157', 'Hunter42', 'escape_time')

% pb.plot_spikes_per_event('PV106', 'Hunter4', 'escape_time', ymax=20)
pb.plot_avg_event_spikes_per_cluster('PV106', 'Hunter15', 'escape_time')

%% errors

% ERROR: PV106,Hunter5; Output argument "camTrigs" (and possibly others) not assigned a value in the execution with "OERecording/getCamerasTrigger" function.
% ERROR: PV106,Hunter19; The entered channel number does not exist in the recording!
% ERROR: PV106,Hunter20; The entered channel number does not exist in the recording!
% ERROR: PV106,Hunter21; The entered channel number does not exist in the recording!
% ERROR: PV106,Hunter34; The entered channel number does not exist in the recording!

%% convert to binary for spike sorting

animal_id = 'PV106';
pb = predictionBreak;
recNames = pb.get_recs_for_animal(animal_id);
for i=1:numel(recNames)
    pb.sa.setCurrentRecording(sprintf('Animal=PV106,recNames=%s', recNames{i}));
    OE = pb.sa.currentDataObj;
    outputName = [OE.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'];
    OE.convert2Binary(outputName);
end

%%
foo_regev()
%%
fclose('all')
EscT = readtable('/media/sil3/Data/regev/circle_escape_trials.csv');
SA = sleepAnalysis('/media/sil3/Data/Regev/brainStatesWake.xlsx');
%%

chosenAnimal = 'PV161';
% plotAllSpectrograms(SA, EscT, chosenAnimal, 7)
recNames = unique(EscT.rec_name(strcmp(EscT.animal_id,chosenAnimal)));
runSingleCircleAnalysis(SA, EscT, chosenAnimal, 'Hunter4', 7, false)

% Hunter5 - bad OE recording
% Hunter10 - bad oe event times (3 seconds before)
% Hunter15 - oscillation effects
% Hunter16 - bad oe event times (0.1 sec before)


%%

EscT = readtable('/media/sil3/Data/regev/circle_escape_trials.csv');
G = findgroups(EscT.animal_id, EscT.rec_name);
subTables = splitapply(@(rows){EscT(rows,:)}, (1:height(EscT))', G);
res = cell(5, height(EscT));
k = 1;
for i = 1:numel(subTables)
    s = subTables{i};
    animalId = s{1, "animal_id"}{1};
    recName = s{1, "rec_name"}{1};
    if animalId ~= "PV106" || recName ~= "Hunter15"
        continue
    end
    SA = sleepAnalysis('/media/sil3/Data/Regev/brainStatesWake.xlsx');
    SA.setCurrentRecording(sprintf('Animal=%s,recNames=%s', animalId, recName));
    SA.getFilters();
    arena = SA.getArenaCSVs(1);
    
    t_activ_screen = getScreenDiodeActivations(arena, SA.currentDataObj, trial_id+1);
    


    for j=1:height(s)
        trial_id = s{j, "trial_id"};
        dt = datetime(s{j, 'escape_time'}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'Asia/Jerusalem');
        escape_timestamp = posixtime(dt);
        event_frame = SA.getVideoFrames(arena.videoFrames, escape_timestamp);
        event_frame = event_frame - arena.IRFrames(1);
        event_oe_time = arena.oeCamTrigs(event_frame);
        try
            
    
            delta = (t_activ_screen(end) - event_oe_time) / 1000;
            dt = datetime(s{j, 'exit_time'}, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS', 'TimeZone', 'Asia/Jerusalem');
            exit_timestamp = posixtime(dt);
            delta_exp = exit_timestamp - escape_timestamp;
            exit_frame = SA.getVideoFrames(arena.videoFrames, exit_timestamp);
            exit_frame = exit_frame - arena.IRFrames(1);
            exit_oe_time = arena.oeCamTrigs(exit_frame);
    
            [V, ~] = SA.currentDataObj.getData(7, event_oe_time-1000, 2000);
            HPV = SA.filt.FH.getFilteredData(V);
            res{1, k} = animalId;
            res{2, k} = recName;
            res{3, k} = trial_id;
            res{4, k} = squeeze(V);
            res{5, k} = squeeze(HPV);
            k = k + 1;
        catch
            k = k + 1;
        end
    end
end
emptyCols = all(cellfun(@isempty, res), 1);
res = res(:, ~emptyCols);

%%

T = 0.05:0.05:2000;
cols = 2;
rows = ceil(width(res)/cols);
figure;
tiledlayout(rows, cols, "TileSpacing","tight","Padding","tight");
for i=1:width(res)
    nexttile();
    plot(T, res{5, i});
    xline(1000)
    % ylim([-0.1, 0.1])
    % ylim([-0.01, 0.01])
    set(gca, "xtick", [])
end


%%

animal = "PV106";
SA = sleepAnalysis('/media/sil3/Data/Regev/brainStatesWake.xlsx');
animalRecs = findCircleFlipRecs(SA, animal);
flipsT = table();
for i = 1:numel(animalRecs)
    SA.setCurrentRecording(sprintf('Animal=%s,recNames=%s', animal, animalRecs{i}));
    arena = SA.getArenaCSVs(1);
    newRow = table(animal, {animalRecs{i}}, height(arena.appEvents), 'VariableNames', {'animal', 'rec', 'n_flips'});
    flipsT = [flipsT; newRow];  % Append the new row
    % fprintf('%s,%s: %d', animal, animalRecs{i}, height(arena.appEvents));
end

%%
T = arena.bugs;
% Step 1: Sort table by time (optional but recommended)
T = sortrows(T, 'Timestamps');
% Step 2: Compute time differences between consecutive rows
dt = diff(T.Timestamps);  % convert to seconds
% Step 3: Identify start of new group where difference > 1 second
isNewGroup = [true; dt > 1];
% Step 4: Assign group numbers
groupID = cumsum(isNewGroup);
% Step 5: Optionally, add group info to table
T.Group = groupID;
% Step 6: Split into cell array of tables (one per group)
groups = splitapply(@(rows){T(rows,:)}, (1:height(T))', groupID);

cols = 5;
rows = ceil(numel(groups) / cols);
for i = 1:numel(groups)
    G = groups{i};  % Extract the i-th group (table)
    radius = arena.trialsData.circle_radius(i);
    x = (G.x - arena.trialsData.circle_position{i}(1)) / radius;
    y = (G.y - arena.trialsData.circle_position{i}(2)) / radius;
    subplot(rows,cols,i)
    plot(x, y)
    xlim([-1.2, 1.2])
    ylim([-1.2, 1.2])
    set(gca, 'YDir', 'reverse')
end


%%
T = readtable('/media/sil3/Data/Pogona_Vitticeps/PV106/2_Hunter/PV106_Trial3_Hunter1_CF_2025-07-23_11-41-49/block1/app_events.csv', 'Delimiter', ',');


%%
plot_id = 1;
plotWin = OE.samplingFrequency(1)*0.2;
for i=1:height(arena.startTrigSh)
    tStart = arena.startTrigSh(i);
    tEnd = arena.endTrigSh(i);
    [V_diode, t_diode] = OE.getAnalogData(1, tStart-300, tEnd-tStart+300);
    t_diode = t_diode + tStart-300;
    V_diode = squeeze(V_diode);
    
    below = V_diode < 3.8e5;
    d = diff([0; below; 0]);  % pad with 0s at edges
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    lengths = ends - starts + 1;
    validStarts = starts(lengths >= 100);
    if i < 5
        subplot(4,1,plot_id)
        [V, t] = OE.getData(8, t_diode(validStarts(1))-2000, 4000);
        plot(t, squeeze(V))
        ylim([-0.5, 0.5])
        xline(2000, '--r')
        plot_id = plot_id + 1;
    end
    % plot(t_diode, V_diode)
    % for j=1:length(validStarts)
    %     xline(t_diode(validStarts(j)),'--r')
    % end
    % break;
end

%%
[V, t] = OE.getAnalogData(1, 0, 100000);
V = squeeze(V);
plot(t,V);


%%
% OE = OERecording('/media/sil3/Data/Pogona_Vitticeps/PV157/1_Hunter/PV157_Trial7_Hunter5_2024-03-07_14-12-21/Record Node 101');
OE = OERecording('/media/sil3/Data/Pogona_Vitticeps/PV157/1_Hunter/PV157_Trial44_Hunter20_2024-03-17_13-38-48/Record Node 101/');

%%
data = readtable('/media/sil3/Data/Regev/circle_ends.csv');

%%
data_ = data(~strcmp(data.rec, 'Hunter59'),:);
% data_ = data(strcmp(data.description, 'following')&~strcmp(data.rec, 'Hunter59'),:);
SA = sleepAnalysis('/media/sil3/Data/Regev/brainStatesWake.xlsx');
V = zeros(1,height(data_), 2000*20);
for j=1:height(data_)
    SA.setCurrentRecording(sprintf('Animal=PV157,recNames=%s', data_{j,'rec'}{1}));
    arena = SA.getArenaCSVs(1);
    event_frame = SA.getVideoFrames(arena.videoFrames, data_{j,'time'});
    event_frame = event_frame - arena.IRFrames(1);
    event_oe_time = arena.oeCamTrigs(event_frame);

    [V(1,j,:), t] = SA.currentDataObj.getData(18, event_oe_time-1000, 2000);
end

%%
SA.setCurrentRecording(sprintf('Animal=PV157,recNames=%s', 'Hunter86'));
arena = SA.getArenaCSVs(1);
event_frame = SA.getVideoFrames(arena.videoFrames, 1712757661.38900);
event_frame = event_frame - arena.IRFrames(1);
event_oe_time = arena.oeCamTrigs(event_frame)

%% 
% add to table the bug_type, bug_speed and number of circle trial so far
% extract for each of the trials the bug appearance and exit time
% compare the signal in high frequency for each of the 3 eventa
% look for signal patterns that appear after these events
% maybe use unsupervised temporal clustering or select manually

%%
figure
idx = strcmp(data_.description, 'following');
% imagesc(squeeze(V(1,idx,:)));
imagesc(squeeze(V));
colorbar()

%% avergae plot
figure
plot(squeeze(mean(V, 2)))

%%
plot_spectral(V)
%%
plot_spectral(V(:,strcmp(data_.description, 'following'),:))
%%
figure;
tiledlayout(ceil(height(data_)/2), 2, 'TileSpacing', 'tight', 'Padding', 'tight')
nexttile
plot(t, squeeze(V))
xline(t(length(t)/2), '--r')
title(sprintf('%s %d', data_{j,'rec'}{1}, data_{j,'time'}))
ylim([-500, 500])
if j < height(data_) - 1 
    xticks('')
end
%%

F=filterData(20000);
F.highPassCutoff=200;
F.lowPassCutoff=2500;
F=F.designBandPass;
VF=F.getFilteredData(V);
%%
figure
plotShifted(squeeze(VF)','verticalShift', 130);


%%
T = OE.getTrigger;

%% 

timeSeriesViewer(SA.currentDataObj)
%%
t0 = 473844.05;
[V_uV,t_ms]=OE.getData(8, t0, 25000);
V_uV=squeeze(V_uV)';
camT = OE.getCamerasTrigger(7);
camT = camT - t0;


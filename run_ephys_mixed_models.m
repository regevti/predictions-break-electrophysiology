function results = run_ephys_mixed_models(options)
%RUN_EPHYS_MIXED_MODELS Rebuild Figure 15 spike/LFP statistics.
%   This analysis is deliberately separate from figure6.m. It writes the
%   exact long-form tables passed to each model, model summaries, fitted
%   model objects, and a good-only reversal sensitivity analysis.

    arguments
        options.cacheDir char = ''
        options.outputDir char = ''
        options.runGamma logical = true
        options.runGoodOnlySensitivity logical = true
    end

    scriptDir = fileparts(mfilename('fullpath'));
    projectDir = fileparts(scriptDir);
    if isempty(options.cacheDir)
        options.cacheDir = fullfile(projectDir, 'cache');
    end
    if isempty(options.outputDir)
        options.outputDir = fullfile(projectDir, 'plots', 'mixed_models');
    end
    if ~isfolder(options.outputDir)
        mkdir(options.outputDir);
    end

    specs = condition_specs();
    results = struct();
    for i = 1:numel(specs)
        spec = specs(i);
        fprintf('\nBuilding %s spike model table...\n', spec.name);
        spikeTable = build_spike_table(options.cacheDir, spec);
        tableFile = fullfile(options.outputDir, sprintf('spikes_%s_model_table.csv', spec.name));
        writetable(spikeTable, tableFile);

        fprintf('Fitting %s Poisson GLME (%d rows)...\n', spec.name, height(spikeTable));
        [mdl, summary, fitLog] = fit_spike_model(spikeTable, spec.name);
        writetable(summary, fullfile(options.outputDir, sprintf('spikes_%s_model_summary.csv', spec.name)));
        writelines(fitLog, fullfile(options.outputDir, sprintf('spikes_%s_model.txt', spec.name)));
        save(fullfile(options.outputDir, sprintf('spikes_%s_model.mat', spec.name)), ...
             'mdl', 'summary', '-v7.3');

        results.spikes.(spec.name).table = spikeTable;
        results.spikes.(spec.name).model = mdl;
        results.spikes.(spec.name).summary = summary;

        if strcmp(spec.name, 'reverse_direction') && options.runGoodOnlySensitivity
            sensitivity = reverse_unit_sensitivity(spikeTable);
            writetable(sensitivity, fullfile(options.outputDir, ...
                       'spikes_reverse_direction_unit_sensitivity.csv'));
            results.spikes.(spec.name).unitSensitivity = sensitivity;

            goodTable = spikeTable(spikeTable.cluster_group=='good',:);
            fprintf('Fitting reverse_direction good-only Poisson GLME (%d rows)...\n', height(goodTable));
            [goodMdl, goodSummary, goodFitLog] = fit_spike_model(goodTable, 'reverse_direction_good_only');
            writetable(goodSummary, fullfile(options.outputDir, ...
                       'spikes_reverse_direction_good_only_model_summary.csv'));
            writelines(goodFitLog, fullfile(options.outputDir, ...
                       'spikes_reverse_direction_good_only_model.txt'));
            save(fullfile(options.outputDir, 'spikes_reverse_direction_good_only_model.mat'), ...
                 'goodMdl', 'goodSummary', '-v7.3');
            results.spikes.(spec.name).goodOnlyTable = goodTable;
            results.spikes.(spec.name).goodOnlyModel = goodMdl;
            results.spikes.(spec.name).goodOnlySummary = goodSummary;
        end
    end

    if options.runGamma
        for i = 1:numel(specs)
            spec = specs(i);
            fprintf('\nBuilding %s gamma model table...\n', spec.name);
            gammaTable = build_gamma_table(options.cacheDir, spec);
            tableFile = fullfile(options.outputDir, sprintf('gamma_%s_model_table.csv', spec.name));
            writetable(gammaTable, tableFile);

            fprintf('Fitting %s Gaussian LME (%d rows)...\n', spec.name, height(gammaTable));
            [mdl, summary, fitLog] = fit_gamma_model(gammaTable, spec.name);
            writetable(summary, fullfile(options.outputDir, sprintf('gamma_%s_model_summary.csv', spec.name)));
            writelines(fitLog, fullfile(options.outputDir, sprintf('gamma_%s_model.txt', spec.name)));
            save(fullfile(options.outputDir, sprintf('gamma_%s_model.mat', spec.name)), ...
                 'mdl', 'summary', '-v7.3');

            results.gamma.(spec.name).table = gammaTable;
            results.gamma.(spec.name).model = mdl;
            results.gamma.(spec.name).summary = summary;
        end
    end

    save(fullfile(options.outputDir, 'ephys_mixed_model_results.mat'), 'results', '-v7.3');
end


function specs = condition_specs()
    specs(1).name = 'circular_to_linear';
    specs(1).eventName = 'escape_time';
    specs(1).animals = {'PV106', 'PV157', 'PV153', 'PV126'};
    specs(1).spikeRecs = { ...
        {'Hunter4','Hunter7','Hunter15','Hunter17'}, ...
        {'Hunter35','Hunter43','Hunter51','Hunter54','Hunter87'}, ...
        {'Hunter12','Hunter13','Hunter16','Hunter17','Hunter18','Hunter19', ...
         'Hunter20','Hunter24','Hunter28','Hunter43','Hunter45','Hunter59','Hunter60'}, ...
        {'Hunter10','Hunter11','Hunter12','Hunter22','Hunter45','Hunter50','Hunter51'}};

    specs(2).name = 'reverse_direction';
    specs(2).eventName = 'flip_time';
    specs(2).animals = {'PV106', 'PV143', 'PV153'};
    specs(2).spikeRecs = { ...
        {'Hunter17','Hunter32'}, ...
        {'Hunter20','Hunter21'}, ...
        {'Hunter10','Hunter16','Hunter60'}};
end


function T = build_spike_table(cacheDir, spec)
    animal = strings(0,1);
    session = strings(0,1);
    unit = strings(0,1);
    eventIndex = zeros(0,1);
    epoch = strings(0,1);
    spikeCount = zeros(0,1);
    windowDuration = zeros(0,1);
    clusterGroup = strings(0,1);

    for ai = 1:numel(spec.animals)
        animalId = spec.animals{ai};
        recs = spec.spikeRecs{ai};
        for ri = 1:numel(recs)
            recName = recs{ri};
            cacheFile = fullfile(cacheDir, sprintf( ...
                'spikes_%s_%s_%s_1_1_1_50_1.mat', animalId, recName, spec.eventName));
            assert(isfile(cacheFile), 'Missing spike cache: %s', cacheFile);
            D = load(cacheFile, 'res');
            res = D.res;

            for ui = 1:width(res)
                rate = double(res{2,ui});
                t = double(res{3,ui});
                t = t(:)';
                assert(~isempty(rate) && size(rate,2) == numel(t), ...
                       'Malformed rate matrix in %s, unit column %d.', cacheFile, ui);
                binWidth = median(diff(t));
                [preMask, postMask, preDuration, postDuration] = prediction_break_epoch_masks(t);
                assert(sum(preMask) == sum(postMask), ...
                       'Unequal spike epoch bin counts in %s.', cacheFile);
                assert(abs(sum(preMask)*binWidth - preDuration) < 1e-9 && ...
                       abs(sum(postMask)*binWidth - postDuration) < 1e-9, ...
                       'Spike bins do not span the requested one-second epochs in %s.', cacheFile);

                preCount = sum(rate(:,preMask),2) * binWidth;
                postCount = sum(rate(:,postMask),2) * binWidth;
                assert(all(isfinite(preCount) & isfinite(postCount)), ...
                       'Non-finite spike counts in %s, unit column %d.', cacheFile, ui);
                assert(all(abs(preCount-round(preCount)) < 1e-8) && ...
                       all(abs(postCount-round(postCount)) < 1e-8), ...
                       'Rate-to-count conversion was not integral in %s.', cacheFile);
                preCount = round(preCount);
                postCount = round(postCount);

                nTrials = size(rate,1);
                rows = 2*nTrials;
                clusterId = double(res{1,ui});
                group = string(res{6,ui});
                animal(end+1:end+rows,1) = animalId;
                session(end+1:end+rows,1) = recName;
                unit(end+1:end+rows,1) = string(clusterId);
                eventIndex(end+1:end+rows,1) = repelem((1:nTrials)',2);
                epoch(end+1:end+rows,1) = repmat(["pre";"post"],nTrials,1);
                spikeCount(end+1:end+rows,1) = reshape([preCount postCount].',[],1);
                windowDuration(end+1:end+rows,1) = repmat([preDuration;postDuration],nTrials,1);
                clusterGroup(end+1:end+rows,1) = group;
            end
        end
    end

    T = table(categorical(animal), categorical(session), categorical(unit), ...
              eventIndex, categorical(epoch,["pre","post"]), spikeCount, ...
              windowDuration, log(windowDuration), categorical(clusterGroup), ...
              'VariableNames', {'animal','session','unit','event_index','epoch', ...
                                'spike_count','window_duration','log_duration','cluster_group'});
end


function [mdl, summary, fitLog] = fit_spike_model(T, condition)
    formulas = { ...
        'spike_count ~ epoch + (1|animal) + (1|animal:session) + (1|animal:session:unit)', ...
        'spike_count ~ epoch + (1|animal:session) + (1|animal:session:unit)'};
    dropped = {'none', 'animal random intercept; full model failed'};
    [mdl, formulaUsed, droppedTerm, fitWarning] = fit_with_fallback( ...
        @(formula) fitglme(T, formula, 'Distribution','Poisson', 'Link','log', ...
                           'Offset',T.log_duration, 'FitMethod','Laplace'), ...
        formulas, dropped);

    [randomLabels, randomVariances] = random_variances(mdl);
    mu = fitted(mdl);
    poissonObservationVariance = log(1 + 1/max(mean(mu),eps));
    pearsonResidual = residuals(mdl, 'ResidualType','Pearson');
    overdispersionRatio = sum(pearsonResidual.^2) / mdl.DFE;
    summary = model_summary(mdl, condition, formulaUsed, droppedTerm, fitWarning, ...
        'Poisson', 'log', randomLabels, randomVariances, poissonObservationVariance, T);
    summary.overdispersion_ratio = overdispersionRatio;
    fitLog = compose_model_log(mdl, summary, ...
        "Poisson observation-level variance uses log(1 + 1/mean fitted count) on the latent scale. " + ...
        "overdispersion_ratio is Pearson chi-square divided by residual degrees of freedom.");
end


function T = build_gamma_table(cacheDir, spec)
    animal = strings(0,1);
    session = strings(0,1);
    trial = strings(0,1);
    trialId = zeros(0,1);
    epoch = strings(0,1);
    gammaDb = zeros(0,1);

    for ai = 1:numel(spec.animals)
        animalId = spec.animals{ai};
        cacheFile = fullfile(cacheDir, sprintf( ...
            'event_signals_%s_circle_%s_NaN_1_1_1_1_1.mat', animalId, spec.eventName));
        assert(isfile(cacheFile), 'Missing gamma cache: %s', cacheFile);
        D = load(cacheFile, 'res');
        res = D.res;
        for ri = 1:width(res)
            recName = char(res{1,ri});
            V = double(res{2,ri});
            fs = double(res{3,ri});
            S = event_aligned_spectrogram(V, fs, ...
                'FreqRange',[1 150], 'WindowSec',0.2, 'StepSec',0.01, ...
                'PreWindow',[-1 0], 'FullWindow',[-1 1], 'Plot',false);
            band = squeeze(mean(S.P_event_dB(S.f>=40 & S.f<=100,:,:),1));
            if isvector(band)
                band = band(:);
            end
            [preMask, postMask] = prediction_break_epoch_masks(S.t);
            pre = mean(band(preMask,:),1,'omitnan')';
            post = mean(band(postMask,:),1,'omitnan')';
            assert(all(isfinite(pre) & isfinite(post)), ...
                   'Non-finite gamma epoch means in %s, %s.', cacheFile, recName);

            nTrials = numel(pre);
            ids = double(res{5,ri});
            ids = ids(:);
            assert(numel(ids) == nTrials, 'Trial-ID mismatch in %s, %s.', cacheFile, recName);
            rows = 2*nTrials;
            animal(end+1:end+rows,1) = animalId;
            session(end+1:end+rows,1) = recName;
            trial(end+1:end+rows,1) = string(repelem(ids,2));
            trialId(end+1:end+rows,1) = repelem(ids,2);
            epoch(end+1:end+rows,1) = repmat(["pre";"post"],nTrials,1);
            gammaDb(end+1:end+rows,1) = reshape([pre post].',[],1);
        end
    end

    T = table(categorical(animal), categorical(session), categorical(trial), trialId, ...
              categorical(epoch,["pre","post"]), gammaDb, ...
              'VariableNames', {'animal','session','trial','trial_id','epoch','gamma_db'});
end


function [mdl, summary, fitLog] = fit_gamma_model(T, condition)
    formulas = { ...
        'gamma_db ~ epoch + (1|animal) + (1|animal:session) + (1|animal:session:trial)', ...
        'gamma_db ~ epoch + (1|animal:session) + (1|animal:session:trial)'};
    dropped = {'none', 'animal random intercept; full model failed'};
    [mdl, formulaUsed, droppedTerm, fitWarning] = fit_with_fallback( ...
        @(formula) fitlme(T, formula, 'FitMethod','REML'), formulas, dropped);

    [randomLabels, randomVariances, residualVariance] = random_variances(mdl);
    summary = model_summary(mdl, condition, formulaUsed, droppedTerm, fitWarning, ...
        'Gaussian', 'identity', randomLabels, randomVariances, residualVariance, T);
    summary.overdispersion_ratio = NaN;
    fitLog = compose_model_log(mdl, summary, ...
        "Gamma model includes a trial random intercept to preserve pre/post pairing.");
end


function [mdl, formulaUsed, droppedTerm, fitWarning] = fit_with_fallback(fitter, formulas, dropped)
    errors = strings(0,1);
    for i = 1:numel(formulas)
        lastwarn('');
        try
            mdl = fitter(formulas{i});
            [warningMessage, warningId] = lastwarn;
            formulaUsed = formulas{i};
            droppedTerm = dropped{i};
            if isempty(warningMessage)
                fitWarning = "none";
            else
                fitWarning = string(warningId) + ": " + string(warningMessage);
            end
            return
        catch ME
            errors(end+1,1) = string(formulas{i}) + " -> " + string(ME.identifier) + ": " + string(ME.message);
        end
    end
    error('predictionBreak:MixedModelFailure', ...
          'All mixed-model structures failed:\n%s', strjoin(errors,newline));
end


function [labels, variances, residualVariance] = random_variances(mdl)
    [psi, residualVariance] = covarianceParameters(mdl);
    groups = mdl.Formula.GroupingVariableNames;
    labels = strings(numel(groups),1);
    variances = zeros(numel(groups),1);
    for i = 1:numel(groups)
        labels(i) = strjoin(string(groups{i}), ':');
        value = double(psi{i});
        variances(i) = value(1,1);
    end
end


function S = model_summary(mdl, condition, formulaUsed, droppedTerm, fitWarning, ...
                           distribution, link, randomLabels, randomVariances, residualVariance, T)
    C = mdl.Coefficients;
    coefIdx = strcmp(C.Name, 'epoch_post');
    assert(sum(coefIdx)==1, 'Could not identify the epoch_post coefficient.');
    beta = C.Estimate(coefIdx);
    lower = C.Lower(coefIdx);
    upper = C.Upper(coefIdx);

    variance = containers.Map(cellstr(randomLabels), num2cell(randomVariances));
    animalVar = map_value(variance, 'animal');
    sessionVar = map_value(variance, 'animal:session');
    if strcmp(distribution,'Poisson')
        lowestVar = map_value(variance, 'animal:session:unit');
    else
        lowestVar = map_value(variance, 'animal:session:trial');
    end
    totalVar = animalVar + sessionVar + lowestVar + residualVariance;
    boundaryLabels = randomLabels(randomVariances < 1e-8);
    if isempty(boundaryLabels)
        boundaryTerms = "none";
    else
        boundaryTerms = strjoin(boundaryLabels, ';');
    end

    nUnits = NaN;
    if ismember('unit',T.Properties.VariableNames)
        nUnits = height(unique(T(:,{'animal','session','unit'}),'rows'));
    end
    if ismember('event_index',T.Properties.VariableNames)
        nTrials = height(unique(T(:,{'animal','session','event_index'}),'rows'));
    else
        nTrials = height(unique(T(:,{'animal','session','trial'}),'rows'));
    end

    if strcmp(link,'log')
        ratio = exp(beta);
        ratioLower = exp(lower);
        ratioUpper = exp(upper);
    else
        ratio = NaN;
        ratioLower = NaN;
        ratioUpper = NaN;
    end

    S = table(string(condition), string(formulaUsed), string(distribution), string(link), ...
        string(mdl.FitMethod), string(droppedTerm), string(fitWarning), boundaryTerms, ...
        height(T), numel(categories(T.animal)), ...
        height(unique(T(:,{'animal','session'}),'rows')), nUnits, nTrials, ...
        beta, C.SE(coefIdx), C.tStat(coefIdx), C.DF(coefIdx), C.pValue(coefIdx), lower, upper, ...
        ratio, ratioLower, ratioUpper, ...
        animalVar, sessionVar, lowestVar, residualVariance, ...
        animalVar/totalVar, sessionVar/totalVar, lowestVar/totalVar, ...
        animalVar/totalVar, (animalVar+sessionVar)/totalVar, ...
        (animalVar+sessionVar+lowestVar)/totalVar, ...
        'VariableNames', {'condition','formula','distribution','link','fit_method', ...
        'dropped_term','fit_warning','boundary_terms','n_rows','n_animals','n_sessions','n_units','n_trials', ...
        'epoch_estimate','epoch_se','epoch_statistic','epoch_df','epoch_p','epoch_lower','epoch_upper', ...
        'epoch_ratio','epoch_ratio_lower','epoch_ratio_upper', ...
        'animal_variance','session_variance','lowest_level_variance','residual_variance', ...
        'animal_vpc','session_vpc','lowest_level_vpc', ...
        'icc_same_animal','icc_same_session','icc_same_lowest_level'});
end


function value = map_value(M, key)
    if isKey(M,key)
        value = M(key);
    else
        value = 0;
    end
end


function lines = compose_model_log(mdl, summary, note)
    lines = [ ...
        "MODEL OUTPUT";
        string(evalc('disp(mdl)'));
        "SUMMARY";
        string(evalc('disp(summary)'));
        "NOTE";
        note;
        "Convergence status is recorded as the MATLAB fit warning. R2025a does not expose a public ConvergenceInfo property for these model classes." ...
    ];
end


function result = reverse_unit_sensitivity(T)
    labels = ["all_non_noise","good_only"];
    result = table();
    for li = 1:numel(labels)
        if labels(li) == "good_only"
            X = T(T.cluster_group=='good',:);
        else
            X = T;
        end
        [G,animal,session,unit,epoch] = findgroups(X.animal,X.session,X.unit,X.epoch);
        counts = splitapply(@sum,X.spike_count,G);
        exposure = splitapply(@sum,X.window_duration,G);
        U = table(animal,session,unit,epoch,counts./exposure, ...
                  'VariableNames',{'animal','session','unit','epoch','rate'});
        W = unstack(U,'rate','epoch');
        animals = categories(W.animal);
        for ai = 1:numel(animals)
            A = W(W.animal==animals{ai},:);
            paired = isfinite(A.pre) & isfinite(A.post);
            pre = A.pre(paired);
            post = A.post(paired);
            [~,p,ci,stats] = ttest(post,pre,'Tail','right');
            try
                pSignrank = signrank(post,pre,'Tail','right');
            catch
                pSignrank = NaN;
            end
            row = table(labels(li), string(animals{ai}), numel(pre), mean(post-pre), ...
                        stats.tstat, stats.df, p, ci(1), ci(2), pSignrank, ...
                        'VariableNames', {'cluster_set','animal','n_units','mean_change_hz', ...
                        't_statistic','df','p_right_ttest','ci_lower','ci_upper','p_right_signrank'});
            result = [result; row]; %#ok<AGROW>
        end
    end
end

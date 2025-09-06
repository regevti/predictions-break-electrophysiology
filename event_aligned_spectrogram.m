function S = event_aligned_spectrogram(V, fs, varargin)
    % EVENT_ALIGNED_SPECTROGRAM  Event-locked LFP/MUA time–frequency map (−1 to +1 s)
    %
    %   S = event_aligned_spectrogram(V, fs, Name, Value, ...)
    %
    %   Inputs
    %     V        : n_events x n_samples matrix. Rows are trials/events,
    %                columns are samples spanning [pre, post] around t=0.
    %     fs       : sampling rate (Hz)
    %
    %   Name–Value options (all optional)
    %     'WindowSec'   : Spectrogram window length in seconds (default 0.2)
    %     'StepSec'     : Step between windows in seconds (default 0.02)
    %     'FreqRange'   : [fmin fmax] Hz to display (default [1 150])
    %     'PreWindow'   : [t0 t1] seconds for baseline (default [-1 0])
    %     'FullWindow'  : [t0 t1] seconds covered by V (default [-1 1])
    %     'NFFT'        : FFT length (default: nextpow2(windowSamples))
    %     'Plot'        : true/false to make figure (default true)
    %
    %   Output struct S
    %     .P_event      : power per event (freq x time x events), linear units
    %     .P_event_dB   : baseline-normalized power per event (dB)
    %     .P_mean_dB    : mean across events (freq x time), dB
    %     .t            : time vector (s, 0 = event)
    %     .f            : frequency vector (Hz)
    %     .baselineIdx  : logical index of baseline time bins
    %
    %   Example
    %     % V: n_events x 40000, fs = 30000 (covers [-1, +1] s)
    %     S = event_aligned_spectrogram(V, 30000, 'FreqRange', [2 120], ...
    %                                   'WindowSec', 0.2, 'StepSec', 0.01);
    
    % -------------------- Parse inputs --------------------
    p = inputParser;
    addParameter(p, 'WindowSec', 0.2, @(x)isscalar(x)&&x>0);
    addParameter(p, 'StepSec',   0.02, @(x)isscalar(x)&&x>0);
    addParameter(p, 'FreqRange', [1 150], @(x)isvector(x)&&numel(x)==2);
    addParameter(p, 'PreWindow', [-1 0], @(x)isvector(x)&&numel(x)==2);
    addParameter(p, 'FullWindow',[-1 1], @(x)isvector(x)&&numel(x)==2);
    addParameter(p, 'NFFT',      [], @(x)isempty(x)||isscalar(x));
    addParameter(p, 'Plot',      true, @(x)islogical(x)||ismember(x,[0 1]));
    addParameter(p, 'Title', 'Event-aligned spectrogram (baseline-normalized, dB)', @(x)ischar(x))
    addParameter(p, 'Ax',      false, @(x)islogical(x)||ismember(x,[0 1]));
    addParameter(p, 'Clim', [], @(x)isvector(x)||isempty(x));
    parse(p, varargin{:});
    opt = p.Results;
    
    [nEvents, nSamples] = size(V);
    
    % Time axis for V (assumed uniform over FullWindow)
    t = linspace(opt.FullWindow(1), opt.FullWindow(2), nSamples);
    
    % Spectrogram params
    wSamp   = max(4, round(opt.WindowSec * fs));        % at least a few samples
    stepS   = max(1, round(opt.StepSec   * fs));
    overlap = max(0, wSamp - stepS);
    if isempty(opt.NFFT), nfft = 2^nextpow2(wSamp); else, nfft = opt.NFFT; end
    
    % Preallocate (first pass to get f, tt)
    [xS, f, tt] = spectrogram(V(1,:), hamming(wSamp,'periodic'), overlap, nfft, fs);
    P1 = abs(xS).^2;                 % power
    % Align spectrogram time to event time (centered w.r.t. original t window)
    % spectrogram returns times as centers relative to start of the signal.
    % Map to absolute event-locked time:
    t0 = t(1); Tspan = t(end)-t(1);
    tt_event = t0 + (tt / (t(end)-t(1))) * Tspan;  %#ok<*NASGU>  (kept for clarity)
    % Simpler & accurate mapping:
    tt_event = linspace(opt.FullWindow(1), opt.FullWindow(2), numel(tt));
    
    P_event   = zeros(numel(f), numel(tt), nEvents, 'like', P1);
    
    % -------------------- Compute per-event spectrograms --------------------
    for ev = 1:nEvents
        x = V(ev,:);
        [X, ~, ~] = spectrogram(x, hamming(wSamp,'periodic'), overlap, nfft, fs);
        P_event(:,:,ev) = abs(X).^2;
    end
    
    % Restrict freq range for display
    fmask = f >= opt.FreqRange(1) & f <= opt.FreqRange(2);
    f = f(fmask);
    P_event = P_event(fmask,:,:);
    
    % Baseline normalization (dB relative to mean pre)
    baselineIdx = tt_event >= opt.PreWindow(1) & tt_event <= opt.PreWindow(2);
    P_event_dB = zeros(size(P_event), 'like', P_event);
    for ev = 1:nEvents
        Pb = mean(P_event(:,baselineIdx,ev), 2);      % mean over baseline time
        Pb(Pb<=0) = eps;                               % guard
        P_event_dB(:,:,ev) = 10*log10(bsxfun(@rdivide, P_event(:,:,ev), Pb));
    end
    P_mean_dB = mean(P_event_dB, 3);
    
    % -------------------- Package output --------------------
    S = struct();
    S.P_event    = P_event;
    S.P_event_dB = P_event_dB;
    S.P_mean_dB  = P_mean_dB;
    S.t          = tt_event;
    S.f          = f;
    S.baselineIdx= baselineIdx;
    
    % -------------------- Plot --------------------
    if opt.Plot
        if ~opt.Ax
            figure();
        end
        clim = opt.Clim;
        if isempty(clim)
            clim = [min(S.P_mean_dB(:)) max(S.P_mean_dB(:))];
        end
        imagesc(S.t, S.f, S.P_mean_dB, clim);
        set(gca,'YDir','normal');
        xlabel('Time (s)'); ylabel('Frequency (Hz)');
        title(opt.Title,'Interpreter','none');
        colormap(parula); cb = colorbar; ylabel(cb,'dB vs. baseline');
        hold on; xline(0,'k-','LineWidth',1.2);
        % Shade baseline window
        yl = ylim;
        patch([opt.PreWindow(1) opt.PreWindow(2) opt.PreWindow(2) opt.PreWindow(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
              'FaceAlpha',0.25, 'EdgeColor','none');
        uistack(findobj(gca,'Type','Line','-and','Color',[0 0 0]), 'top');
    end
end
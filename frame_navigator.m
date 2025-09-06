function frame_navigator(video_path, start_frame, event_frames, event_labels, title)
% FRAME_NAVIGATOR  Simple GUI to browse video frames starting at a given frame.
%
%   frame_navigator(video_path, start_frame, event_frames, event_labels)
%     - video_path    : path to video file
%     - start_frame   : (default 1) starting frame index
%     - event_frames  : (optional) vector of frame numbers to jump to via dropdown
%     - event_labels  : (optional) labels for events (cellstr/string), same length as event_frames
%     - title         : (optional) title for the figure
%
%   Keys: ←/→ step 1, PgUp/PgDn step 10, Home/End to ends.

    arguments
        video_path (1,1) string
        start_frame (1,1) double {mustBeInteger,mustBePositive} = 1
        event_frames double {mustBeInteger,mustBePositive} = []
        event_labels = []   % accept string array, cellstr, or empty; validate below
        title char = 'Frame Navigator'
    end

    v = VideoReader(video_path);
    estTotalFrames = max(1, floor(v.Duration * v.FrameRate));   % estimate

    % Sanitize events now that we know total frames
    if ~isempty(event_frames)
        event_frames = unique(event_frames(:).');
        event_frames = event_frames(event_frames>=1 & event_frames<=estTotalFrames);
    end

    % ----- Event labels handling -----
    if isempty(event_frames)
        labels = {};
    else
        % Normalize event_labels to a cellstr if provided; otherwise auto-generate
        if nargin < 4 || isempty(event_labels)
            labels = arrayfun(@(k) sprintf('Event %d',k), 1:numel(event_frames), 'UniformOutput', false);
        else
            % Convert to cellstr (column) and validate length
            if isstring(event_labels)
                labels = cellstr(event_labels(:));
            elseif iscell(event_labels)
                labels = cellfun(@(x) string(x), event_labels(:), 'UniformOutput', false);
                labels = cellstr(string(labels)); % ensure cellstr
            else
                error('event_labels must be a string array, cell array of char, or empty.');
            end
            if numel(labels) ~= numel(event_frames)
                error('event_labels must have the same length as event_frames.');
            end
        end
    end

    % Clamp start
    start_frame = min(max(start_frame, 1), estTotalFrames);

    % --- UI ---
    fig = uifigure('Name',title,'Position',[789 983 1920 1080]);
    gl  = uigridlayout(fig,[3 8]);
    gl.RowHeight   = {'1x', 35, 35};
    gl.ColumnWidth = {'1x','fit','fit','fit','fit','1x','fit','fit'};

    ax  = uiaxes(gl); ax.Layout.Row = 1; ax.Layout.Column = [1 8];
    ax.Toolbar.Visible = 'off'; ax.Interactions = []; axis(ax,'off');

    btnPrev = uibutton(gl,'Text','⟵ Prev','ButtonPushedFcn',@(s,e) step(-1));
    btnPrev.Layout.Row = 2; btnPrev.Layout.Column = 2;

    btnNext = uibutton(gl,'Text','Next ⟶','ButtonPushedFcn',@(s,e) step(+1));
    btnNext.Layout.Row = 2; btnNext.Layout.Column = 3;

    lblF   = uilabel(gl,'Text','Frame:'); lblF.Layout.Row=2; lblF.Layout.Column=4;
    editF  = uieditfield(gl,'numeric','Limits',[1 estTotalFrames], ...
                         'Value',start_frame,'RoundFractionalValues',true, ...
                         'ValueDisplayFormat','%.0f');
    editF.Layout.Row=2; editF.Layout.Column=5;
    btnGo  = uibutton(gl,'Text','Go','ButtonPushedFcn',@(s,e) gotoFrame(editF.Value));
    btnGo.Layout.Row=2; btnGo.Layout.Column=6;

    % Events dropdown
    lblEvt = uilabel(gl,'Text','Event:');
    lblEvt.Layout.Row = 2; lblEvt.Layout.Column = 7;

    ddEvt = uidropdown(gl);
    ddEvt.Layout.Row = 2; ddEvt.Layout.Column = 8;

    if isempty(event_frames)
        ddEvt.Items = {"No events"};
        ddEvt.ItemsData = nan;
        ddEvt.Value = ddEvt.ItemsData(1);
        ddEvt.Enable = 'off';
    else
        ddEvt.Items = labels(:)';        % your provided (or default) labels
        ddEvt.ItemsData = event_frames;  % frames as the data
        % Default to first event or nearest ≥ start_frame
        ix = find(event_frames >= start_frame, 1, 'first');
        if isempty(ix), ix = 1; end
        ddEvt.Value = event_frames(ix);
        ddEvt.ValueChangedFcn = @(s,e) gotoFrame(e.Value);
    end

    sld = uislider(gl,'Limits',[1 estTotalFrames],'Value',start_frame, ...
                   'MajorTicks',[],'MinorTicks',[]);
    sld.Layout.Row = 3; sld.Layout.Column = [1 8];
    sld.ValueChangingFcn = @(s,e) gotoFrame(round(e.Value));
    sld.ValueChangedFcn  = @(s,e) gotoFrame(round(e.Value));

    lblInfo = uilabel(gl,'Text',''); lblInfo.Layout.Row=2; lblInfo.Layout.Column=1;

    % Enable arrow keys
    fig.KeyPressFcn = @(~,ev) onKey(ev);

    % state
    curFrame = start_frame;
    show(curFrame);

    % ------- nested helpers -------
    function onKey(ev)
        switch ev.Key
            case {'leftarrow'} , step(-1);
            case {'rightarrow'}, step(+1);
            case {'pageup'}    , step(-10);
            case {'pagedown'}  , step(+10);
            case {'home'}      , gotoFrame(1);
            case {'end'}       , gotoFrame(estTotalFrames);
        end
    end

    function step(delta)
        gotoFrame(curFrame + delta);
    end

    function gotoFrame(f)
        f = min(max( round(f), 1), estTotalFrames);
        if f ~= curFrame
            curFrame = f;
            show(curFrame);
        else
            syncUI();
        end
    end

    function show(f)
        v.CurrentTime = (f-1) / v.FrameRate;
        try
            img = readFrame(v);
        catch
            v.CurrentTime = max(0, (f-1)/v.FrameRate - 0.5/v.FrameRate);
            if hasFrame(v), img = readFrame(v); else, return; end
        end
        imshow(img,'Parent',ax,'Border','tight');
        axis(ax,'image'); axis(ax,'off');
        syncUI();
    end

    function syncUI()
        sld.Value   = curFrame;
        editF.Value = curFrame;
        curTime = (curFrame-1)/v.FrameRate;
        lblInfo.Text = sprintf('%.3f s  |  Frame %d / %d  |  %.2f fps', ...
                               curTime, curFrame, estTotalFrames, v.FrameRate);

        % Keep dropdown in sync if current frame is exactly an event
        if ~isempty(event_frames) && any(ddEvt.ItemsData == curFrame)
            ddEvt.Value = curFrame;
        end
        drawnow limitrate
    end
end

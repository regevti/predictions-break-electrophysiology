function play_video_segment(videoPath, startFrame, numFrames)
% play_video_segment Play a segment of a video
%
%   play_video_segment(videoPath, startFrame, numFrames)
%   - videoPath:   path to the video file (string/char)
%   - startFrame:  first frame index to play (integer, >=1)
%   - numFrames:   number of frames to play (integer)

    % Open video
    v = VideoReader(videoPath);

    % Validate start frame
    if startFrame < 1 || startFrame > floor(v.Duration * v.FrameRate)
        error('Start frame is out of bounds.');
    end

    % Calculate start time
    startTime = (startFrame-1) / v.FrameRate;
    v.CurrentTime = startTime;

    % Create figure
    hFig = figure('Name','Video Segment Player','NumberTitle','off');

    % Loop through frames
    for k = 1:numFrames
        if hasFrame(v)
            frame = readFrame(v);
            imshow(frame, 'Border', 'tight');
            title(sprintf('Frame %d', startFrame + k - 1));
            pause(1/v.FrameRate); % match video speed
        else
            break;
        end

        % Stop if figure closed
        if ~isvalid(hFig)
            break;
        end
    end
end
function vid_file = find_video_files(path_in, cam_name)
    % Get parent directory
    parent_dir = fileparts(path_in);
    % Recursive search for files starting with "right"
    files_avi = dir(fullfile(parent_dir, '**', sprintf('%s*.avi',cam_name)));
    files_mp4 = dir(fullfile(parent_dir, '**', sprintf('%s*.mp4',cam_name)));
    files = [files_avi; files_mp4];
    % Filter out directories (keep only files)
    files = files(~[files.isdir]);
    % Display found file paths
    vid_file = nan;
    if isempty(files)
        fprintf('No files starting with "%s" found in %s\n', cam_name, parent_dir);
    elseif numel(files) > 1
        fprintf('More than 1 file starting with "%s" found in %s\n', cam_name, parent_dir);
    else
        vid_file = fullfile(files(1).folder, files(1).name);
    end
end

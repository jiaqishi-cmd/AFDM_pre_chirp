function rootDir = find_afdm_root(startDir)
%FIND_AFDM_ROOT Locate the matlab_afdm source root from a script folder.
%   Experiment scripts may live in subfolders. This helper walks upward
%   until it finds setup_paths.m, which marks the matlab_afdm root.

    if nargin < 1 || isempty(startDir)
        startDir = fileparts(mfilename('fullpath'));
    end

    rootDir = startDir;
    while true
        if exist(fullfile(rootDir, 'setup_paths.m'), 'file') == 2
            return;
        end

        parentDir = fileparts(rootDir);
        if strcmp(parentDir, rootDir)
            error('Could not locate matlab_afdm root from %s.', startDir);
        end
        rootDir = parentDir;
    end
end

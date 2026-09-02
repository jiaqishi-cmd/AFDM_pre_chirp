function setup_paths(rootDir)
%SETUP_PATHS Add AFDM source folders to the MATLAB path.
%   setup_paths() resolves folders relative to this file.
%   setup_paths(rootDir) resolves folders relative to rootDir.

    if nargin < 1 || isempty(rootDir)
        rootDir = fileparts(mfilename('fullpath'));
    end

    addpath(rootDir);
    receiveDir = fullfile(rootDir, 'receive');
    if exist(receiveDir, 'dir')
        addpath(receiveDir);
    end

    toolsDir = fullfile(rootDir, 'tools');
    if exist(toolsDir, 'dir')
        addpath(toolsDir);
    end

    experimentsDir = fullfile(rootDir, 'experiments');
    if exist(experimentsDir, 'dir')
        addpath(genpath(experimentsDir));
    end
end

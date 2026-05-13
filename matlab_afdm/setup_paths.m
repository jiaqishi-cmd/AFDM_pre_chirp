function setup_paths(rootDir)
%SETUP_PATHS Add AFDM source folders to the MATLAB path.
%   setup_paths() resolves folders relative to this file.
%   setup_paths(rootDir) resolves folders relative to rootDir.

    if nargin < 1 || isempty(rootDir)
        rootDir = fileparts(mfilename('fullpath'));
    end

    addpath(rootDir);
    experimentsDir = fullfile(rootDir, 'experiments');
    if exist(experimentsDir, 'dir')
        addpath(genpath(experimentsDir));
    end
end

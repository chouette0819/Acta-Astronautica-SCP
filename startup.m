% MATLAB startup script for this project
% Adds the GA-STM folder (if present) to the path recursively.

if exist('GA-STM','dir')
    addpath(genpath('GA-STM'));
end

% You can add more project-specific paths below as needed.


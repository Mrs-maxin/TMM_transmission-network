% TEMPLATE SCAFFOLD: imported from project-template-lite; not recovered research implementation.
% Cloud Python and original project MATLAB files remain pending recovery.
function ensure_dir(dirPath)
%ENSURE_DIR Create a directory if it does not already exist.

dirPath = string(dirPath);
if strlength(dirPath) == 0
    error("TMM_Transmission_utils:io:EmptyDir", "Directory path is empty.");
end

if ~exist(dirPath, "dir")
    mkdir(dirPath);
end
end

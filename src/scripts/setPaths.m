function [dirname, basedir] = setPaths()
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
basedir = getAbsPath(fullfile(dirname, '../..'));
if length(size(basedir)) > 1
    basedir = basedir{1};
end
warning('off', 'MATLAB:dispatcher:nameConflict')
addpath(genpath(fullfile(basedir,'src')));
addpath(genpath(fullfile(basedir,'data')));

function absPath = getAbsPath(obj)
getAbsFolderPath = @(y) string(unique(arrayfun(@(x) x.folder, dir(y), 'UniformOutput', false)));
getAbsFilePath = @(y) string(arrayfun(@(x) fullfile(x.folder, x.name), dir(y), 'UniformOutput', false));
if isfolder(obj)
    absPath = getAbsFolderPath(obj);
elseif isfile(obj)
    absPath = getAbsFilePath(obj);
else
    error('The specified object does not exist.');
end
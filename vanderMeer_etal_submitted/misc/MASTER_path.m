% path used for Alyssa's SWR project
restoredefaultpath;

% get machine name
if ispc
    machinename = getenv('COMPUTERNAME');
    filesep = '\';
else
    machinename = getenv('HOSTNAME');
    filesep = '/';
end

% retrive location of github path
switch machinename
    case 'ISIDRO'
        github_path = 'C:\Users\mvdm\Documents\GitHub';
    case {'EQUINOX','BERGKAMP'}
        github_path = 'D:\My_Documents\GitHub';
end

addpath(genpath(cat(2,github_path,'\vandermeerlab\code-matlab\shared')));
addpath(genpath(cat(2,github_path,'\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze')));
addpath(genpath(cat(2,github_path,'\vandermeerlab\code-matlab\tasks\Replay_Analysis')));

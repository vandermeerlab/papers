function traj_out = GenerateAHVTrajectory(cfg_in)
% function traj_out = GenerateAHVTrajectory(cfg_in)
%
%

cfg_def = [];
cfg_def.mode = 1;
cfg_def.t = [0 5];
cfg_def.dt = 1/60;
cfg_def.ahv = 40; % deg/s

cfg = ProcessConfig(cfg_def,cfg_in);

tvec = cfg.t(1):cfg.dt:cfg.t(2);
nSamples = length(tvec);

switch cfg.mode
    case 0 % constant
    
        ahv = cfg.ahv*ones(size(tvec));
        traj_out = tsd(tvec,ahv);
        
    case 1 % predefined example
        
        split = floor(nSamples/2);
        ahv(1:split) = cfg.ahv; % zig
        ahv(split+1:nSamples) = -cfg.ahv; % zag
        traj_out = tsd(tvec,ahv);
        
    otherwise
        error('Mode not yet implemented.');
end
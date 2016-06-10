function sdf_out = GenerateSDFfromTC(cfg_in,hd,tc)
% function sdf_out = GenerateSDFfromTC(cfg_in,hd,tc)
%

cfg_def = [];
cfg_def.mode = 'poisson'; % {'interp','poisson'}

cfg = ProcessConfig(cfg_def,cfg_in);

nCells = size(tc.tc,1);

for iC = nCells:-1:1
    
    sdf_out(iC,:) = interp1(tc.xbin,tc.tc(iC,:),hd.data,'linear');
    
    % hack to handle values out of range (should do circular interpolation)
    nan_idx = find(isnan(sdf_out(iC,:)));
    sdf_out(iC,nan_idx) = interp1(tc.xbin,tc.tc(iC,:),hd.data(nan_idx),'nearest','extrap');
    
    if strcmp(cfg.mode,'poisson') % generate Poisson rates with lambda taken from TC
    
        len = length(hd.data);
        poiss_rng = rand(1,len);

        sdf_out(iC,:) = poissinv(poiss_rng,sdf_out(iC,:));
    
    end
    
end

sdf_out = tsd(hd.tvec,sdf_out);
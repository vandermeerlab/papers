function [fake_data] = AMPX_Naris_fake_data_power(cfg, nEvts)
%% creates a series of power gradents with some noise and a set of random gradients ('"conrol")
%          Inputs: 
%           - nEvts: N how many events to generate for lg, hg, lg_ran,
%           hg_ran
%          Outputs: 
%           - fake_power: [struct] for each band .power_distrib{length(nEvt)} = 8x8 gradients
%                                                .power_distrib_avg = 8x8 average across events
% EC - 2016-11-27

fake_data = [];

% make a parametric distribution of values with a mean of 1 and an SD of
% 0.1
rng('shuffle')
r = 1+0.1*randn(10000,1);
% hist(r,100)
bands = {'lg', 'hg'};
for iband = 1:2
    for ievt = nEvts:-1:1
        [x y] = meshgrid(1:2);
        grad = ones(2,2);
        grad(2,2) = 10*iband^2;
        %     surf(x,y,grad_100)
        [Xq,Yq] = meshgrid(1:1/8:2);
        Vq = interp2(x,y,grad,Xq,Yq);
        
        Vq = Vq(2:end, 2:end);
        Xq = Xq(2:end, 2:end);
        Yq = Yq(2:end, 2:end);
        %     surf(Xq,Yq,Vq_w_noise)
        % add some parametric noise
        rng('shuffle')
        ran_id = randperm(1000, 64);
        % create the noisey gradient
        if cfg.noise
        fake_data.(bands{iband}).power.power_distrib(:,:,ievt) = Vq .* reshape(r(ran_id), 8,8); 
        else
                    fake_data.(bands{iband}).power.power_distrib(:,:,ievt) = Vq; 

        end
        % create a random "control" gradient
        
        r = 2.5+1*randn(10000,1);
        ran_id = randperm(1000, 64);
        
        fake_data.([bands{iband} '_ran']).power.power_distrib(:,:,ievt) = reshape(r(ran_id), 8,8);
    end
end

%% do the same for some spindles

for ievt = nEvts:-1:1
        [x y] = meshgrid(1:2);
        grad = ones(2,2);
        grad(1,1) = 5;
        %     surf(x,y,grad_100)
        [Xq,Yq] = meshgrid(1:1/8:2);
        Vq = interp2(x,y,grad,Xq,Yq);
        
        Vq = Vq(2:end, 2:end);
        Xq = Xq(2:end, 2:end);
        Yq = Yq(2:end, 2:end);
        %     surf(Xq,Yq,Vq)
        % add some parametric noise
        rng('shuffle')
        ran_id = randperm(1000, 64);
        % create the noisey gradient
        fake_data.spindles.power.power_distrib(:,:,ievt) = Vq .* reshape(r(ran_id), 8,8);    
end

%% get the averages across evts for each band 
bands = {'lg', 'hg', 'lg_ran', 'hg_ran', 'spindles'};
for iband  = 1:length(bands)
    fake_data.(bands{iband}).power.power_distrib_avg = nanmean(fake_data.(bands{iband}).power.power_distrib,3); 
    [m, n, z] = size(fake_data.(bands{iband}).power.power_distrib);
    fake_data.(bands{iband}).power.power_distrib = squeeze(mat2cell(fake_data.(bands{iband}).power.power_distrib, m, n, ones(1,z)))';
    fprintf(['\n ' bands{iband} '.power.power_distrib nEvents = ' num2str(z) '\n'])
%     fake_data.(bands{iband}).power_distrib_avg
end

%% check some numbers

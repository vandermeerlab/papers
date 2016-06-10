function data_out = WillDataLoader(cfg_in,sno)
% function data_out = WillDataLoader(cfg,sno)
%
% data inventory and loader

cfg_def.dt = 1/60;
cfg_def.columns = 'FGHIJ';
cfg_def.fd = 'D:\My_Documents\Dropbox\projects\HDfit\data';
cfg_def.ATIshift = 0; % shift HD data forward relative to spikes by this amount

cfg = ProcessConfig(cfg_def,cfg_in);

p = pwd;
cd(cfg.fd);

switch sno
    case 1 % WB84 4-15
        sess = {'std','laser'};
        fn = {'WB84 4-15 s1 ST7 c1 ST8 c1 std.xls','WB84 4-15 s3 ST7 c1 ST8 c1 darklaseron.xlsx'};
        sn = {'WB84 4-15 s1 ST7 c1 ST8 c1 std.','Sheet1'};
        nCells = 2;
        
    case 2 % WB84 3-23
        sess = {'std','laser'};
        fn = {'WB84 3-23 s1 ST6c1 std.xls','WB84 3-23 s3 ST6c1 darklaseron.xls'};
        sn = {'WB84 3-23 s1 ST6c1 std.txt','WB84 3-23 s3 ST6c1 darklaseron.'};
        nCells = 1;
        
    case 3 % WB84 3-26
        sess = {'std','laser'};
        fn = {'WB84 3-26 s1 ST3c1 ST5c1 ST8c1 std.xls','WB84 3-26 s3 ST3c1 ST5c1 ST8c1 darklaseron.xls'};
        sn = {'WB84 3-26 s1 ST3c1 ST5c1 ST8c1 ','WB84 326 s3 ST3c1 ST5c1 ST8c1 d'};
        nCells = 3;
        
    case 4 % WB89 7-1
        sess = {'std','laser'};
        fn = {'WB89 7-1p s1 ST7c1 std.xls','WB89 7-1p s3 ST7c1 darklaseron.xls'};
        sn = {'WB89 7-1p s1 ST7c1 std.txt','WB89 7-1p s3 ST7c1 darklaseron.'};
        nCells = 1;
        
    case 5 % WB89 7-2
        sess = {'std','laser'};
        fn = {'WB89 7-2 s1 ST7c1 std.xls','WB89 7-2p s2 ST7c1 darklaseron.xls'};
        sn = {'WB89 7-2 s1 ST7c1 std.txt','WB89 7-2p s2 ST7c1 darklaseron.'};
        nCells = 1;
        
    case 6 % WB89 7-7
        sess = {'std','laser'};
        fn = {'WB89 7-7 s1 ST7c1 std.xls','WB89 7-7 s4 ST7c1 darklaseron.xls'};
        sn = {'WB89 7-7 s1 ST7c1 std.txt','WB89 7-7 s4 ST7c1 darklaseron.t'};
        nCells = 1;
        
    case 7 % WB95 8-28
        sess = {'std','laser'};
        fn = {'WB95 8-28p s1 ST7c1 std.xls','WB95 8-28p s3 ST7c1 darklaseron.xls'};
        sn = {'WB95 8-28p s1 ST7c1 std.txt','WB95 8-28p s3 ST7c1 darklaseron'};
        nCells = 1;
        
    case 8 % WB95 9-1
        sess = {'std','laser'};
        fn = {'WB95 9-1 s1 ST1c1 ST8c1 std.xls','WB95 9-1 s3 ST1c1 ST8c1 darklaseron.xls'};
        sn = {'WB95 9-1 s1 ST1c1 ST8c1 std.txt','WB95 9-1 s3 ST1c1 ST8c1 darklas'};
        nCells = 2;
        
    case 9 % WB95 9-7
        sess = {'std','laser'};
        fn = {'WB95 9-6 s1 ST7c1 ST8c2 std.xls','WB95 9-7 s2 ST7c1 ST8c2 darklaseron.xls'};
        sn = {'WB95 9-6 s1 ST7c1 ST8c2 std.txt','WB95 9-7 s2 ST7c1 ST8c2 darklas'};
        nCells = 2;
        
    case 10 % WB85 3-22
        sess = {'std','laser'};
        fn = {'WB85 3-22 s1 ST6c2 std.xls','WB85 3-22 s3 ST6c2 darklaseron.xls'};
        sn = {'WB85 3-22 s1 ST6c2 std.txt','WB85 3-22 s3 ST6c2 darklaseron.'};
        nCells = 1;
        
    case 11 % WB85 3-23
        sess = {'std','laser'};
        fn = {'WB85 3-23 s1 ST8c1c2c3c4 std.xls','WB85 3-23 s3 ST8c1c2c3c4 darklaseron.xls'};
        sn = {'WB85 3-23 s1 ST8c1c2c3c4 std.tx','WB85 3-23 s3 ST8c1c2c3c4 darkla'};
        nCells = 4;
        
    case 12 % WB85 3-26
        sess = {'std','laser'};
        fn = {'WB85 3-26 s1 ST5c1 std.xls','WB85 3-26 s3 ST5c1 darklaseron.xls'};
        sn = {'WB85 3-26 s1 ST5c1 std.txt','WB85 3-26 s3 ST5c1 darklaseron.'};
        nCells = 1;
        
    case 13 % WB76 2-9
        sess = {'std','laser'};
        fn = {'WB76 2-9 s1 ST6c1c2 std.xls','WB76 2-9 s3 ST6c1c2 darklaseron.xls'};
        sn = {'WB76 2-9 s1 ST6c1c2 std.txt','WB76 2-9 s3 ST6c1c2 darklaseron'};
        nCells = 2;
        
    case 14 % WB76 2-9p
        sess = {'std','laser'};
        fn = {'WB76 2-9pm s1 ST6c1c2 std.xls','WB76 2-9pm s3 ST6c1c2 darklaseron.xls'};
        sn = {'WB76 2-9pm s1 ST6c1c2 std.txt','WB76 2-9pm s3 ST6c1c2 darklaser'};
        nCells = 2;
        
    case 15 % WB76 2-10
        sess = {'std','laser'};
        fn = {'WB76 2-10 s3 ST6c1 std2.xls','WB76 2-10 s5 ST6c1 darklaseron.xls'};
        sn = {'WB76 2-10 s3 ST6c1 std2.txt','WB76 2-10 s5 ST6c1 darklaseron.'};
        nCells = 1;
        
end

for iF = 1:length(fn)
   
    data_out.(sess{iF}).fn = fn{iF};
    
    bin_idx = xlsread(fn{iF},sn{iF},'A1:A65536')';
    bin_idx = bin_idx(1:end-cfg.ATIshift);
    
    data_out.(sess{iF}).bin_idx = bin_idx;
     
    for iC = 1:nCells
    
        this_column = cfg.columns(iC);
        this_range = cat(2,this_column,'1:',this_column,'65536');
        data_out.(sess{iF}).obs_fr(iC,:) = xlsread(fn{iF},sn{iF},this_range);
         
    end
    data_out.(sess{iF}).obs_fr = data_out.(sess{iF}).obs_fr(:,1:end-cfg.ATIshift);
       
    
    this_column = cfg.columns(nCells+1);
    this_range = cat(2,this_column,'1:',this_column,'65536');
    obs_hd = xlsread(fn{iF},sn{iF},this_range); obs_hd = obs_hd(1+cfg.ATIshift:end);

    orig_tvec = (bin_idx-1)*cfg.dt;
    data_out.(sess{iF}).obs_hd = tsd(orig_tvec,obs_hd);
    
    %%% HANDLE SOME SPECIAL CASES %%%
    if sno == 3 & iF == 2 % remove some troublesome VT samples which trip up the kalman filters

        data_out.(sess{iF}) = remove_samples(data_out.(sess{iF}),[7878 7879]);
        data_out.(sess{iF}) = remove_samples(data_out.(sess{iF}),[305 306]);
        
    end
    
    if sno == 11 % throw out cells with strange tuning
        data_out.(sess{iF}).obs_fr = data_out.(sess{iF}).obs_fr([1 4],:);
    end
        
end

cd(pwd);


function data = remove_samples(data,samples_to_remove)

data.obs_hd.data(samples_to_remove) = [];
data.obs_hd.tvec(samples_to_remove) = [];
data.obs_fr(:,samples_to_remove) = [];
data.bin_idx(samples_to_remove) = [];
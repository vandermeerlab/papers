%%
cd('D:\My_Documents\Dropbox\projects\Sushant\2014-10-17');
load fd_equinox


%% create
for iFD = 2:length(fd)
    
   cd(fd{iFD});
   
   MakeTrialinfo([]);
   
   %pause;
   
   ExpandTrialinfo([]);
  
    
end

%% collect
for iFD = 1:length(fd)
    
   cd(fd{iFD});
   load(FindFile('*info.mat'));
   
   if iFD == 1
      ALL_trialinfo = trialinfo; 
   else
       fn = fieldnames(trialinfo);
       fprintf('Session %s, nFieldnames = %d\n',fd{iFD},length(fn));
       
       for ifn = 1:length(fn)
          ALL_trialinfo.(fn{ifn}) = cat(2,ALL_trialinfo.(fn{ifn}),trialinfo.(fn{ifn})); 
           
       end
   end
   
   
end

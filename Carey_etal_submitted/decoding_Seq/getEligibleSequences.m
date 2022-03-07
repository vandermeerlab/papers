function seq_out = getEligibleSequences(cfg_in,decSeq)
% seq_out = getEligibleSequences(cfg,decSeq)
%
% Identify and organize set of eligible sequences based on raw detected sequences (output of Generate_DecSeq.m)
%
% INPUTS:
%
% out: variable obtained by loading single-session output file from GenerateDecSeq.m
% 
% OUTPUTS:
%
% seq_out: iv describing eligible sequences, with the following usr fields:
% 
% .side: 1 (left), 2 (right)
% .shuf_count: proportion of identity shuffles also detected as a sequence (averaged across samples)
% .nActiveCells
% .start_loc: integer decoded bin number
% .end_loc: integer decoded bin number
% ...and outputs of linear regression fit (.beta, .pval, .rsq)
%
% CONFIG:
%
% .evt: if specified, use this as candidate SWR events, i.e. only include those candidates that overlap with this input

swapfun = @(x) -x+3; % utility function to access variable from other condition

cfg_select = [];
cfg_select.verbose = 0;

cfg_def = [];
cfg_def.whichSeq = 'all'; % {'all','fwd','bwd','either'}; % which sequences to process?
cfg_def.arms = {'left','right'};
cfg_def.sessions = {'food','water'};
cfg_def.minActiveCells = 4;
cfg_def.minlen = 0.05; % otherwise, minimum length in s
cfg_def.RemoveOverlappingEvents = 1; % remove events that are classified as both L and R -- default is 1
cfg_def.processShuffles = 1;
cfg_def.shuffle_include_p = 0.05; % discard detected interval if more than this fraction of its bins is detected in shuffles

cfg = ProcessConfig(cfg_def,cfg_in);

%% for each arm, restrict sequences by SWR overlap, minimum length, and number of cells active
for iCond = 1:length(cfg.arms)
        
    fprintf('%s, number of raw sequences: %d\n',cfg.arms{iCond},length(decSeq.expCond(iCond).seq_iv.tstart));
    
    nShuffles = length(decSeq.expCond(iCond).seq_iv_shuf);
    
        % only include events that overlap with SWRs if specified
        if isfield(cfg,'evt')
            decSeq.expCond(iCond).seq_iv = IntersectIV([],decSeq.expCond(iCond).seq_iv,cfg.evt);
            
            if cfg.processShuffles
                
                for iShuf = 1:nShuffles
                    decSeq.expCond(iCond).seq_iv_shuf{iShuf} = IntersectIV([],decSeq.expCond(iCond).seq_iv_shuf{iShuf},cfg.evt);
                end
            end
        end
        fprintf('%s, number of sequences after SWR overlap: %d\n',cfg.arms{iCond},length(decSeq.expCond(iCond).seq_iv.tstart));
        
        % only select events with some minimum length
        if ~isempty(cfg.minlen)
            evt_len = decSeq.expCond(iCond).seq_iv.tend-decSeq.expCond(iCond).seq_iv.tstart;
            decSeq.expCond(iCond).seq_iv = SelectIV(cfg_select,decSeq.expCond(iCond).seq_iv,evt_len >= cfg.minlen);
            
            if cfg.processShuffles
                for iShuf = 1:nShuffles
                    evt_len = decSeq.expCond(iCond).seq_iv_shuf{iShuf}.tend - decSeq.expCond(iCond).seq_iv_shuf{iShuf}.tstart;
                    decSeq.expCond(iCond).seq_iv_shuf{iShuf} = SelectIV(cfg_select,decSeq.expCond(iCond).seq_iv_shuf{iShuf},evt_len >= cfg.minlen);
                end
            end
        end
        fprintf('%s, number of sequences after minimum length: %d\n',cfg.arms{iCond},length(decSeq.expCond(iCond).seq_iv.tstart));
        
        % only select events with minimum number of cells
        if ~isempty(decSeq.expCond(iCond).seq_iv.tstart)
            decSeq.expCond(iCond).seq_iv = AddNActiveCellsIV(cfg_select,decSeq.expCond(iCond).seq_iv,decSeq.expCond(iCond).decS);
            keep_idx = decSeq.expCond(iCond).seq_iv.usr.nActiveCells >= cfg.minActiveCells;
            decSeq.expCond(iCond).seq_iv = SelectIV(cfg_select,decSeq.expCond(iCond).seq_iv,keep_idx);
            
            if cfg.processShuffles
                for iShuf = 1:nShuffles
                    if ~isempty(decSeq.expCond(iCond).seq_iv_shuf{iShuf}.tstart)
                        decSeq.expCond(iCond).seq_iv_shuf{iShuf} = AddNActiveCellsIV(cfg_select,decSeq.expCond(iCond).seq_iv_shuf{iShuf},decSeq.expCond(iCond).decS);
                        keep_idx = decSeq.expCond(iCond).seq_iv_shuf{iShuf}.usr.nActiveCells >= cfg.minActiveCells;
                        decSeq.expCond(iCond).seq_iv_shuf{iShuf} = SelectIV(cfg_select,decSeq.expCond(iCond).seq_iv_shuf{iShuf},keep_idx);
                    end
                end
            end
        end
        fprintf('%s, number of sequences after minimum number of cells: %d\n',cfg.arms{iCond},length(decSeq.expCond(iCond).seq_iv.tstart));
        
        
        % remove sequences which don't pass shuffle threshold
        if cfg.processShuffles
            nShuf = length(decSeq.expCond(iCond).seq_iv_shuf);
            shuf_count = decSeq.expCond(iCond).decode_map; % initialize counts to zero
            shuf_count.data = zeros(size(shuf_count.data));
            
            for iShuf = 1:nShuf
                [~,shuf_idx] = restrict(shuf_count,decSeq.expCond(iCond).seq_iv_shuf{iShuf});
                shuf_count.data(shuf_idx) = shuf_count.data(shuf_idx) + 1;
            end
            shuf_count.data = shuf_count.data./nShuf;
            
            add_cfg = []; add_cfg.label = 'shuf_count'; add_cfg.method = 'mean';
            decSeq.expCond(iCond).seq_iv = AddTSDtoIV(add_cfg,decSeq.expCond(iCond).seq_iv,shuf_count);
            
            nRemovedShuf = sum(decSeq.expCond(iCond).seq_iv.usr.shuf_count > cfg.shuffle_include_p); % track number of events removed this way
            decSeq.expCond(iCond).seq_iv = SelectIV(cfg_select,decSeq.expCond(iCond).seq_iv,decSeq.expCond(iCond).seq_iv.usr.shuf_count <= cfg.shuffle_include_p);
        end
        fprintf('%s, number of sequences after shuffle removal: %d\n',cfg.arms{iCond},length(decSeq.expCond(iCond).seq_iv.tstart));
        
        
end % of left/right

%% take stock
this_seq(1) = decSeq.expCond(1);
this_seq(2) = decSeq.expCond(2);

%% remove overlaps if requested
for iCond = 1:length(cfg.arms)
    
    if cfg.RemoveOverlappingEvents
        other_seq = decSeq.expCond(swapfun(iCond)).seq_iv;
        
        this_seq(iCond).seq_iv = DifferenceIV(cfg_select,this_seq(iCond).seq_iv,other_seq);
        % hmm does not handle usr fields?!
        
        if cfg.processShuffles
            for iShuf = 1:nShuffles
                other_seq = decSeq.expCond(swapfun(iCond)).seq_iv_shuf{iShuf};
                this_seq(iCond).seq_iv_shuf{iShuf} = DifferenceIV(cfg_select,this_seq(iCond).seq_iv_shuf{iShuf},other_seq);
            end
        end
    end
    fprintf('%s, number of sequences after overlap removal: %d\n',cfg.arms{iCond},length(this_seq(iCond).seq_iv.tstart));
        
end

%% add fwd/bwd information
for iCond = 1:length(cfg.arms)
    
    [this_seq(iCond).seq_iv.usr.beta,this_seq(iCond).seq_iv.usr.rsq,this_seq(iCond).seq_iv.usr.pval,this_seq(iCond).seq_iv.usr.start_loc,this_seq(iCond).seq_iv.usr.end_loc] = ...
        deal(nan(size(this_seq(iCond).seq_iv.tstart))); % initialize vars to column shapes
    
    this_map = this_seq(iCond).decode_map;
    for iS = 1:length(this_seq(iCond).seq_iv.tstart)
        
        this_iv = restrict(this_map,this_seq(iCond).seq_iv.tstart(iS),this_seq(iCond).seq_iv.tend(iS));
        
        [betas,~,~,~,stat] = regress(this_iv.data',[ones(1,length(this_iv.data)); 1:length(this_iv.data)]');
        
        this_seq(iCond).seq_iv.usr.beta(iS) = betas(2); % regression coefficient
        this_seq(iCond).seq_iv.usr.rsq(iS) = stat(1); % variance explained
        this_seq(iCond).seq_iv.usr.pval(iS) = stat(3); % p-value
        this_seq(iCond).seq_iv.usr.start_loc(iS) = this_iv.data(1);
        this_seq(iCond).seq_iv.usr.end_loc(iS) = this_iv.data(end);
        
    end
 
end

%% now add identity to usr field, and unify
for iCond = 1:length(cfg.arms)
    if ~isempty(this_seq(iCond).seq_iv.usr.nActiveCells)
        this_seq(iCond).seq_iv.usr.side = iCond*ones(size(this_seq(iCond).seq_iv.usr.nActiveCells));
    else
        this_seq(iCond).seq_iv.usr.side = iCond*ones(0,1);
    end
end
seq_out = UnionIV([],this_seq(1).seq_iv,this_seq(2).seq_iv);
function iv_out = InvertIV(iv_in,tsd_in)
% function iv_out = InvertIV(iv_in,tsd_in)
% 
% 

t0 = tsd_in.tvec(1);
t1 = tsd_in.tvec(end);

iv_out = iv;
iv_out.tstart = cat(1,t0,iv_in.tend);
iv_out.tend = cat(1,iv_in.tstart,t1);


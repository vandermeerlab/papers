function sig_out = AUX_shuffle_phases(sig_in)
% function sig_out = AUX_shuffle_phases(sig_in)
%
% use inverse fourier transform to generate signal with randomly permuted
% phases
%
% MvdM 2016-01-11

md = abs(fft(sig_in));
phi = angle(fft(sig_in));

len = floor(length(phi)/2);
phi_half = phi(1:len);
phi_shuf_half = phi_half(randperm(len)); % permute phases

if mod(length(phi),2) % odd length
    phi_shuf = cat(2,phi_shuf_half,0,-fliplr(phi_shuf_half));
else
    phi_shuf = cat(2,phi_shuf_half,-fliplr(phi_shuf_half));
end

% fr = ifft(md.*exp(1i*phi)); % original signal

sig_out = ifft(md.*exp(1i*phi_shuf),'symmetric'); % reconstruct signal with shuffled phases


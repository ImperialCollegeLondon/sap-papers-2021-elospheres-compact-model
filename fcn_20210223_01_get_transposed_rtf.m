function[d_t] = fcn_20210223_01_get_transposed_rtf(ema,src_az,src_inc,nfft)

h = ema.getImpulseResponseForSrc(src_az,src_inc); %[nFIR, nChan]
H = v_rfft(h,nfft,1); 
d = H./H(:,ema.refChan); % [nFreq,nChan]
d_t = d.'; %[nChan,nFreq]
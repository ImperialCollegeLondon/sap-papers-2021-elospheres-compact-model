function[R] = fcn_20170421_02_diffuse_pow_to_psd_matrix(ema,fs,nfft,az,inc,pow,quad)
% given the power arriving from each specified direction (and an optional
% quadrature scalar for each direction) simulate the microphone signals 
% which would be experienced in a diffuse field and from these compute an
% mvdr beamformer


% ema: handle of instantiated ElobesMicArray object
az = az(:);
inc = inc(:);
pow = pow(:);
if nargin>6 && ~isempty(quad)
    quad = quad(:);
else
    quad = ones(size(az));
end

ema.prepareData(fs);

h = ema.getImpulseResponseForSrc(az(:),inc(:)); %[len_filt,nSensors,nDirections]
H = v_rfft(h,nfft,1); % [nFreq,nChan,nDOA]
[nFreq,nChan,nDOA] = size(H);

per_doa_weight = sqrt(pow).*quad; %[nDOA 1]

% 3D [nFreq,nChan,nDOA] -> 2D [nDOA nFreqxnChan]
reshape(permute(H,[3 1 2]),nDOA,[]);
% quadrature weighting
Hw=bsxfun(@times,reshape(permute(H,[3 1 2]),nDOA,[]),per_doa_weight);  % weighted H
% 2D [nDOA nFreqxnChan] -> 3D [nFreq,nChan,nDOA]
Hw=permute(reshape(Hw,nDOA,nFreq,nChan),[2 3 1]);


R = bsxfun(@times,permute(Hw,[2 4 1 3]),conj(permute(Hw,[4 2 1 3]))); %[nChan nChan nFreq nDOA]
R = mean(R,4); %[nChan nChan nFreq]


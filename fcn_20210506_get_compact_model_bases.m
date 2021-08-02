function[out] = fcn_20210506_get_compact_model_bases(in_params)

required_fields = {'fs','stft_params','ema','doas','iso_grid'};
if any(~ismember(required_fields,fieldnames(in_params)))
    error('in_params must contain all required fields:\n%s',...
        sprintf('\t%s\n',string(required_fields.')))
end
required_fields2={'azi_step','inc_step'};
if any(~ismember(required_fields2,fieldnames(in_params.iso_grid)))
    error('in_params.iso_grid must contain all required fields:\n%s',...
        sprintf('\t%s\n',string(required_fields2.')))
end

% doas
h_a = in_params.ema.getImpulseResponseForSrc(in_params.doas(:,1),in_params.doas(:,2));
H_a = rfft(h_a,in_params.stft_params.nfft,1);
[nFreq,nChan,nDOA] = size(H_a);

% isotropic quadrature weights for spherical grid
n_inc = round(180/in_params.iso_grid.inc_step);
n_az = round(360/in_params.iso_grid.azi_step);
az_v = (0:n_az-1) * (2*pi/n_az); %vector of azimuths
inc_v=(1:2:2*n_inc)*pi*0.5/n_inc;%vector of inclinations
vq=(n_inc-abs(n_inc+1-2*(1:n_inc))).^(-1).*exp(-1i*(inc_v+0.5*pi));
weight_inc=(-2*sin(inc_v).*real(fft(vq).*exp(-1i*(0:n_inc-1)*pi/n_inc))/n_inc);
%----
az_frac = 2*pi/n_az;
% get a full grid
[inc_g,az_g] = ndgrid(inc_v(:),az_v(:));
quad_g = repmat(1/(4*pi) * az_frac * weight_inc(:),1,length(az_v));
% sanity check
if sum(quad_g(:))-1 > 1e-14,error('quadrature weights don''t sum to one'),end
% for spherically isotropic noise the power is equal in all directions
pow_g = ones(size(inc_g));


% ---- covariance matrices ----
% covariance for isotropic soundfield (quadrature-weighted grid)
out.Riso = fcn_20170421_02_diffuse_pow_to_psd_matrix(in_params.ema,...
    in_params.stft_params.fs,in_params.stft_params.nfft,az_g,inc_g,pow_g,quad_g); % [nSenors,nSensors,nFreq]

% covariance for each direction
out.Rdoas = bsxfun(@times,permute(H_a,[2 4 1 3]),conj(permute(H_a,[4 2 1 3]))); %[nChan nChan nFreq nDOA]

% covariance assuming same power in each mic
out.Rwhite = repmat(eye(nChan),1,1,nFreq); %[nChan nChan nFreq]
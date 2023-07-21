function ap = tfpower_exp1(EEG,t2)

% downsampled time points
times2save = -3500:25:500;

tidx = dsearchn(EEG.times',times2save');
[~,t2_idx(1)] = min(abs(t2(1)-times2save));
[~,t2_idx(2)] = min(abs(t2(2)-times2save));

% soft-coded parameters
freq_range  = [1 30]; % extract only these frequencies (in Hz)
num_frex    = 30;     % number of frequencies between lowest and highest

% set range for variable number of wavelet cycles
range_cycles = [1 12];

% set up convolution parameters
wavtime = -2:1/EEG.srate:2;
frex    = linspace(freq_range(1),freq_range(2),num_frex);
nKern   = length(wavtime);
nData   = EEG.pnts*EEG.trials;
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;

% number of cycles
numcyc = linspace(range_cycles(1),range_cycles(end),num_frex);

% initialize TF matrix
tf = zeros(EEG.nbchan,num_frex,length(tidx),EEG.trials);

% Fourier spectrum of data
dataX = fft( reshape(EEG.data,EEG.nbchan,[]),nConv ,2);

% create wavelets and do TF decomposition in one loop
for fi=1:num_frex

    % create time-domain wavelet
    twoSsquared = 2 * (numcyc(fi)/(2*pi*frex(fi))) ^ 2;
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / twoSsquared );

    % compute fourier coefficients of wavelet and normalize
    cmwX = fft(cmw,nConv);
    cmwX = cmwX / max(cmwX);

    % now loop over channels
    for chani=1:EEG.nbchan

        % second and third steps of convolution
        as = ifft( dataX(chani,:).*cmwX );

        % cut wavelet back to size of data
        as = as(halfwav+1:end-halfwav);
        as = reshape(as,EEG.pnts,EEG.trials);

        powts = abs(as).^2 ;

        % baseline-normalized power time series
        tf(chani,fi,:,:)  = powts(tidx,:);

    end % end channel loop
end % end frequency loop

ap = pow2db(mean(squeeze(mean(squeeze(mean(tf([1,2,7,8,10,11],8:13,t2_idx(1):t2_idx(2),:)))))));



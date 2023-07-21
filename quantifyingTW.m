function [fw,bw]=quantifyingTW(trial,sr)

         f1 = 8;
         f2 = 13;
         power_fft = abs(fftshift(fft2(trial)));
        
        % parameters
        samplingRate = sr;
        numberElectrodes = size(trial,1);
        yElec = 1:numberElectrodes;
        durationSignal = size(trial,2)/samplingRate; 
        dF = 1/durationSignal;
        fx = -samplingRate/2:dF:samplingRate/2-dF;
        fy = (size(trial,1)/2)*linspace(-1,1,size(trial,1));
               
        % calculate BW and FW
        [aBW,bBW] = max(power_fft(fy>0, fx>=f1 & fx<=f2));
        [bwValue,bBW2] = max(aBW);
        bwTempFreq = fx(bBW2+ceil(length(fx)/2)+round(f1/dF));
        bwSpatFreq = fy(bBW(bBW2)+ceil(length(yElec)/2));
        
        [aFW,bFW] = max(power_fft(fy<0, fx>=f1 & fx<=f2));
        [fwValue,bFW2] = max(aFW);
        fwTempFreq = fx(bFW2+ceil(length(fx)/2)+round(f1/dF));
        fwSpatFreq = fy(bFW(bFW2));
        
        % calculate BWss and FWss
        for shuffle_index = 1:100

            trial_shuffle = trial(randperm(size(trial,1)),:);
           
            power_fft_shuffle = abs(fftshift(fft2(trial_shuffle)));

            [aBW,bBW] = max(power_fft_shuffle(fy>0, fx>=f1 & fx<=f2));
            [bwValue_s(shuffle_index),bBW2] = max(aBW);

            [aFW,bFW] = max(power_fft_shuffle(fy<0, fx>=f1 & fx<=f2));
            [fwValue_s(shuffle_index),bFW2] = max(aFW);

        end

        bwValue_shuffle = mean(bwValue_s);
        fwValue_shuffle = mean(fwValue_s);

        bw = 10*log10(bwValue/bwValue_shuffle);
        fw = 10*log10(fwValue/fwValue_shuffle);

end
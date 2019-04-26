function ECoGdata_perblock_LNfilteredDetrended = ...
    sfa_expt4_ECoGrawdata_DetrendLineNoiseFilter(i_block, sub, subs_PreProcSettings, filterOrder, data_ECoGraw_blocks)

    bandstop_freqs = subs_PreProcSettings.(sub).bandStopFreqs; %line noise freqs
    
    fsample = data_ECoGraw_blocks.fsample; %Sample freq in Hz

    ECoGdata_perblock_LNfilteredDetrended = detrend(data_ECoGraw_blocks.trial{i_block}')'; %detrend

        %band reject/stop filter
        for i_bp = 1:length(bandstop_freqs)
            bandStop = [bandstop_freqs(i_bp)-subs_PreProcSettings.(sub).bandStopWidth bandstop_freqs(i_bp)+subs_PreProcSettings.(sub).bandStopWidth];
            [npb,npa] = butter(filterOrder,(bandStop*2)/fsample,'stop');
            ECoGdata_perblock_LNfilteredDetrended = filtfilt(npb,npa,ECoGdata_perblock_LNfilteredDetrended')';    
        end
        
        %highpass filter
        highPassFreq = subs_PreProcSettings.(sub).HighPassFreq; %then high-pass
        [lpb,lpa] = butter(filterOrder,(highPassFreq*2)/fsample,'high');
        ECoGdata_perblock_LNfilteredDetrended = filtfilt(lpb,lpa,ECoGdata_perblock_LNfilteredDetrended')';

%         figure;
%         pwelch(ECoGdata_perblock_LNfilteredDetrended(1,:),[],[],[],512); %filtered data

end

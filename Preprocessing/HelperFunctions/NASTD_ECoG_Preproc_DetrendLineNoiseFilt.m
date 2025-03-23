function ECoGdata_perblock_LNfilteredDetrended = ...
    NASTD_ECoG_Preproc_DetrendLineNoiseFilt...
    (i_block, sub, subs_PreProcSettings, ...
    filterOrder, bandstop_freqs, highPassFreq, ...
    data_ECoGraw_blocks)

%Aim: Detrend, Demean, and Filter block-wise iEEG raw data
fsample = data_ECoGraw_blocks.fsample; %Sample freq in Hz

ECoGdata_perblock_LNfilteredDetrended = ...
    detrend(data_ECoGraw_blocks.trial{i_block}')'; %detrend & demean

%     %visually check
%     i_chan = randi(length(data_ECoGraw_blocks.label),1);
%     figure; hold on
%     plot(data_ECoGraw_blocks.trial{i_block}(i_chan,:),'b-')
%     plot(ECoGdata_perblock_LNfilteredDetrended(i_chan,:),'r--')

%Band-stop filter for line noise
for i_bp = 1:length(bandstop_freqs)
    bandStop = [bandstop_freqs(i_bp)-subs_PreProcSettings.(sub).bandStopWidth ...
        bandstop_freqs(i_bp)+subs_PreProcSettings.(sub).bandStopWidth];
    [npb,npa] = butter(filterOrder,(bandStop*2)/fsample,'stop');
    ECoGdata_perblock_LNfilteredDetrended = filtfilt(npb,npa,ECoGdata_perblock_LNfilteredDetrended')';
end

%Highpass filter
[lpb,lpa] = butter(filterOrder,(highPassFreq*2)/fsample,'high');
ECoGdata_perblock_LNfilteredDetrended = filtfilt(lpb,lpa,ECoGdata_perblock_LNfilteredDetrended')';

%         figure;
%         pwelch(ECoGdata_perblock_LNfilteredDetrended(1,:),[],[],[],512); %filtered data

end

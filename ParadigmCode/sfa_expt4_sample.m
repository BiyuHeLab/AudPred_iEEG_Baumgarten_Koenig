function sfa_expt2_sample(tr_str)

% randomly select an unused series
load series_selection_scale2Hz

tr_ind = tr_str + 1;
rand_ind = ceil(rand*97);

ser     = series_u_Hz{tr_ind}(rand_ind, :);
log_ser = series_u_logHz{tr_ind}(rand_ind, :);

figure; hold on;
plot(log_ser, 'bo-');
ylabel('pitch','FontSize',20);
xlabel('time','FontSize',20);
title({'currently playing tone sequence', ['trend strength = ' num2str(tr_str)]},'FontSize',20);
set(gca,'XTick',[]);
set(gca,'YTick',[]);

drawnow

% initialize psychportaudio
AssertOpenGL;

lowLatency = 1;
InitializePsychSound( lowLatency );

deviceID       = [];
mode           = 1;  % sound playback only
timingPriority = 2;  % take aggressive control of audio device
freq           = []; % allow PsychPortAudio to set the determined optimal value, depending on other settings 
nChannels      = 2;
bufferSize     = 128*2^3;
pahandle       = PsychPortAudio('Open', deviceID, mode, timingPriority, freq, nChannels, bufferSize);

s    = PsychPortAudio('GetStatus', pahandle);
freq = s.SampleRate;

% make sound wave
param.f_sample_inHz  = freq;
param.toneDur_inSecs = .3;
w = series2soundwave(ser, param.toneDur_inSecs, param.f_sample_inHz);

w2 = [w; w];
PsychPortAudio('FillBuffer', pahandle, w2);

mins = min(log_ser);
maxs = max(log_ser);

WaitSecs(.3);
PsychPortAudio('Start', pahandle);

for k = 1:33
t0 = GetSecs;
clf; hold on;
plot(log_ser, 'bo-');
plot([k k], [mins maxs], 'r-');
ylabel('pitch','FontSize',20);
xlabel('time','FontSize',20);
title({'currently playing tone sequence', ['trend strength = ' num2str(tr_str)]},'FontSize',20);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
drawnow;
WaitSecs(.3 - (GetSecs - t0));
end

% WaitSecs(10.3);
PsychPortAudio('Close', pahandle);

end
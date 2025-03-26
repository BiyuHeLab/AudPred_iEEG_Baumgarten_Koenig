function test_tone()

%% set up psychportaudio

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

s    = PsychPortAudio('GetStatus', pahandle)
freq = s.SampleRate

dt = 1/freq;

%% make tone

toneDur_inSecs = 2;
toneFreq = 440;

nSamples = toneDur_inSecs * freq;

ts = 1 : dt : toneDur_inSecs;
w  = sin( 2*pi*toneFreq * ts);

w2 = [w; w];
PsychPortAudio('FillBuffer', pahandle, w2);

PsychPortAudio('Start', pahandle);
WaitSecs(toneDur_inSecs);

PsychPortAudio('Close', pahandle);

end


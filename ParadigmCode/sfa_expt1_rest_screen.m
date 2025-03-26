clear

addpath([pwd '/functions']);
[window keyboardNumbe1r] = openScreen;

%test_trigger(255);


%% get experiment parameters

param = sfa_expt3_getParams(window, [], 0, 1);

Screen('FillRect',  window, param.BGcolor);
Screen('TextSize',  window, param.textSize);
Screen('TextColor', window, param.fontColor);
Screen('TextFont',  window, param.font);

Screen('FillOval', window, param.fontColor, param.fixationRect);
Screen('Flip', window);

%recordValidKeys(GetSecs, Inf, param.keyboardNumber, []);

recordValidKeys(GetSecs, Inf, param.keyboardNumber, param.r1_validKeys);

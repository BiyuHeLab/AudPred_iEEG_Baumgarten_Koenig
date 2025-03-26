function [t_on t_off] = test_trigger(trigger)

config_io;

param.io_address = 16376;
param.io_addressC = 16376;
% param.io_address = hex2dec('3FFA');
% param.io_address = hex2dec('3FF8');
param.io_address = hex2dec('3000');


param.triggerDur_inSecs = .02;

t_on = GetSecs;
% outp(param.io_addressC, 4);
% 
% WaitSecs(param.triggerDur_inSecs);


param.triggerDur_inSecs = .02;

t_on = GetSecs;
outp(param.io_address, trigger);

WaitSecs(param.triggerDur_inSecs);

t_off = GetSecs;
outp(param.io_address, 0)

end
function w = series2soundwave(toneSeries, toneDur_inSecs, f_sample_inHz)


%% sampling frequency
if ~exist('f_sample_inHz','var') || isempty(f_sample_inHz)
    f_sample_inHz = 10000; %8192; % sampling frequency; matlab default is 8000 for sampling sound
end
dt = 1/f_sample_inHz;


%% tone properties
if ~exist('toneDur_inSecs','var') || isempty(toneDur_inSecs)
    toneDur_inSecs = .3;
end
nSamplesInTone = toneDur_inSecs * f_sample_inHz;

nTonesInSeries   = length(toneSeries);
seriesDur_inSecs = toneDur_inSecs * nTonesInSeries;


%% convert the fGn/fBm tone series to a sound wave
w       = [];
phi     = 0; %initial phase set to be 0
phis    = 0;
breakpt = [];
flag    = [];
disable = 0;
debug   = 0;
for f = toneSeries % from the 1st sinnoid to connect all the following sinnoids
    flag(end+1)=0;
    tt  = 0 : dt : toneDur_inSecs;
    aa  = cos(2*pi*f * tt + phi);
    
    if ~disable
        p1  = phi;
        phi = acos(aa(end)); % continuous choice of the phi based on the previous tone
        if sign(sin(2*pi*f * tt(end) + p1)) ~= sign(sin(phi)) % 1st order time derivative continuous
            phi = -phi; 
            flag(end)=1;
        end
        phis(end+1) = phi;
    end
    
    % debug
    if debug
        num_p = 100;
        if ~isempty(w)
            plot([w(end-num_p:end) aa(2:num_p)],'-')  %check the continuity
            hold on
            plot(num_p+1, w(end),'r*')
            title(['flag=' num2str(flag(end))]);

            keyboard
            clf
        end
    end
    
    w = [w aa(2:end)];

%     breakpt(end+1) = length(w)+1;
%     if disable
%         phi=0;
%     end
end

end
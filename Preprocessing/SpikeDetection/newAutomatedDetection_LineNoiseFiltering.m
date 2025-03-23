function d=newAutomatedDetection_LineNoiseFiltering(d,fs,hum_fs)

if nargin<3
    hum_fs=50;
end

if min(size(d))==1
   d=d(:); 
end

R = 1; r = 0.985;


f0 = hum_fs:hum_fs:fs/2; % Hz


for i=1:length(f0)
    b = [1 -2*R*cos(2*pi*f0(i)/fs) R*R];
    a = [1 -2*r*cos(2*pi*f0(i)/fs) r*r];
    for ch=1:size(d,2)
        d(:,ch)=filtfilt(b,a,d(:,ch));
    end
end
end
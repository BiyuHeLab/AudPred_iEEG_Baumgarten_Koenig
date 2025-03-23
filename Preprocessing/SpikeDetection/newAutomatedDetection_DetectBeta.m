
function M=newAutomatedDetection_DetectBeta(data,fs,settings)


M=[];
winsize=settings.muscleArtifactBand_windowSizeSeconds*fs;
noverlap=round(0.5*winsize);

index=1:winsize-noverlap:size(data,1)-winsize+1;
if isempty(index)
    index=1;
    winsize=size(data,2);
end

[bb,aa]=butter(4,2*30/fs);
if license('test','distrib_computing_toolbox')==1
    
    parfor i_chan=1:size(data,2)
        MM=[];
        %     winsize=5*fs;
        %     noverlap=2.5*fs;
        
        for i=1:length(index)
            seg=data(index(i):index(i)+winsize-1,i_chan);
            seg=filtfilt(bb,aa,seg);
            
            a=lpc(seg-mean(seg),settings.muscleArtifactBand_autoregressiveSamples);
            [h,f]=freqz(1,a,[],fs);
            h=abs(h);
            poz=f(diff([0; sign(diff([0; h]))])<0);
            
            MM(i)=sum(poz<25 & poz>settings.muscleArtifactBand_frequency)>0;
        end
        MM=[MM MM(end)];
        M(:,i_chan)=interp1([index size(data,1)],MM,1:size(data,1),'nearest');
    end
    M=logical(M);
else
    
    
    for i_chan=1:size(data,2)
        MM=[];
        %     winsize=5*fs;
        %     noverlap=2.5*fs;
        
        for i=1:length(index)
            seg=data(index(i):index(i)+winsize-1,i_chan);
            seg=filtfilt(bb,aa,seg);
            
            a=lpc(seg-mean(seg),settings.muscleArtifactBand_autoregressiveSamples);
            [h,f]=freqz(1,a,[],fs);
            h=abs(h);
            poz=f(diff([0; sign(diff([0; h]))])<0);
            
            MM(i)=sum(poz<25 & poz>settings.muscleArtifactBand_frequency)>0;
        end
        MM=[MM MM(end)];
        M(:,i_chan)=interp1([index size(data,1)],MM,1:size(data,1),'nearest');
    end
    M=logical(M);
end
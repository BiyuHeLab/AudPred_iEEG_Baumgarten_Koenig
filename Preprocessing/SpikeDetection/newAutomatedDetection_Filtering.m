function df=newAutomatedDetection_Filtering(bandwidth,d,fs,type)

% bandpass filtering
switch type
    case 1
        % IIR-cheby2
        % low pass
        Wp = 2*bandwidth(2)/fs; Ws = 2*bandwidth(2)/fs+ 0.1;
        Rp = 6; Rs = 60;
        [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
        [bl,al] = cheby2(n,Rs,Ws);
        
        % high pass
        Wp = 2*bandwidth(1)/fs; Ws = 2*bandwidth(1)/fs- 0.05;
        Rp = 6; Rs = 60;
        [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
        [bh,ah] = cheby2(n,Rs,Ws,'high');
    case 2
        % IIR butterworth
        [bl,al]=butter(5,2*bandwidth(2)/fs);
        [bh,ah]=butter(5,2*bandwidth(1)/fs,'high');
        
    case 3
        % FIR
        % low pass
        bl = fir1(fs/2,2*bandwidth(2)/fs); al=1;
        
        % high pass
        bh = fir1(fs/2,2*bandwidth(1)/fs,'high'); ah=1;
end


if license('test','distrib_computing_toolbox')==1
    parfor ch=1:size(d,2); df(:,ch)=filtfilt(bh,ah,d(:,ch)); end
    if bandwidth(2)==fs/2; return; end
    parfor ch=1:size(d,2); df(:,ch)=filtfilt(bl,al,df(:,ch)); end
    
else % If you do not have "Parallel Computing Toolbox", only standard for-cycle will be performed
    for ch=1:size(d,2); df(:,ch)=filtfilt(bh,ah,d(:,ch)); end
    if bandwidth(2)==fs/2; return; end
    for ch=1:size(d,2); df(:,ch)=filtfilt(bl,al,df(:,ch)); end
    
end
end


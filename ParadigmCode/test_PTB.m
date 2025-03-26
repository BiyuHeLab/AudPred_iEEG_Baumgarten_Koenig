imgfile = 'example_image1.png';
ima = imread(imgfile, 'png');

sx = 'center';
sy = 'center';
vSpacing = 1.4;

r1_text = 'TESTTESTTEST';
fc = [0,0,0];

try 
    commandwindow;
    nr = max(Screen('Screens'));
    [w, screenRect] = Screen('OpenWindow',nr,0,[],32,2);

    Screen('PutImage', w, ima);
    %Screen('Flip',w);
    DrawFormattedText(w, r1_text, sx, sy, fc, [], [], [], vSpacing);
    Screen('Flip',w);
    while KbCheck; end
    while ~KbCheck; end
    Screen('CloseAll');
    k
catch
    Screen('CloseAll');
    rethrow(lasterror);
end
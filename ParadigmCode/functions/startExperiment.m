function startup = startExperiment(experimentName, isPractice)

disp('working from directory')
disp([pwd ' ...'])
disp(' ')

if IsWin %Psychtoolbox funtion: check if windows or mac version in order to determine thow paths are specified
% if ispc %alternative non-psychtoolbox command
    dataDir = [pwd '\data\'];
else
    dataDir = [pwd '/data/']; %Will be the case, cause mac
end

startupTime = clock;

if isPractice
    currentFile = [experimentName ' PRACTICE' datestr(startupTime,' yyyy-mm-dd HH-MM-SS') '.mat'];
else
    currentFile = [experimentName datestr(startupTime,' yyyy-mm-dd HH-MM-SS') '.mat'];
end
dataFile    = [dataDir currentFile];

startup.dataDir        = dataDir;
startup.dataFile       = dataFile;
startup.startupTime    = startupTime;
startup.startupTimeStr = datestr(startupTime);
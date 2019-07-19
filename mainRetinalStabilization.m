clear all;
close all;
clc;

global TRIALINFO;
global CAMERA;
global STARFIELD;
global STARDATA;
global SCREEN

%% path name and file name
subjectName = '';
fileName = ['retinalStabilization_' subjectName '_' datestr(now,'yymmddHHMM')];
saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

% set keyboard
KbName('UnifyKeyNames'); 
skipKey = KbName('space');
escape = KbName('f1');
leftArror = KbName('LeftArrow');
rightArror = KbName('RightArrow');
upArror = KbName('UpArrow');
cKey = KbName('c');

testMode = 0; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC

%% parameters
coordinateMuilty = 1; % convert cm to coordinate system for moving distance etc.
SCREEN.distance = 60*coordinateMuilty;% cm

if testMode
    SCREEN.widthCM = 34.5*coordinateMuilty; % cm, need to measure in your own PC
    SCREEN.heightCM = 19.7*coordinateMuilty; % cm, need to measure in your own PC
else
    SCREEN.widthCM = 120*coordinateMuilty; % cm
    SCREEN.heightCM = 65*coordinateMuilty; % cm
end

TRIALINFO.repetition = 1;
TRIALINFO.motionType = [4]; % 1: fixation; 2: normal pursuit; 3: simulated pursuit; 4:stabilized pursuit
TRIALINFO.headingDegree = 0 ; % degree, currently unused
TRIALINFO.headingSpeed = 1*coordinateMuilty; % cm/s
TRIALINFO.coherence = 100;
TRIALINFO.fixationSizeD = 0.25;  % degree

% fixation period (fixation point only)          --> 
% pause period (3D cloud appear, not moving)	--> 
% pre movement duration (heading only)          -->
% movement duration (heading with pursuit/fixation)
TRIALINFO.fixationPeriod = 0.5; % second
TRIALINFO.pausePeriod = 0.18; % second
TRIALINFO.preMoveDuration = 0.4; % second
TRIALINFO.moveDuration = 10; % second

TRIALINFO.fixationWindow = 10; % degree
TRIALINFO.pursuitWindow = 30; % degree

TRIALINFO.intertrialInterval = 1; % second

% for motion type 3
TRIALINFO.rotationDegree = [-10 10]; % ¡ã£¬the degree of the star' rotation
TRIALINFO.rotationSpeed = 10;  % ¡ã/s

% for motion type 2 and 4
TRIALINFO.fixMoveDirection = [1 3]; % only for motion type 2 and 4. 1: left to right; 2: constant at the center; 3: right to left;
TRIALINFO.fixationDegree = 4; % degree ¡À to center
TRIALINFO.fixationInitialDegree = 5; % degree ¡À to center
TRIALINFO.fixSpeed = TRIALINFO.rotationSpeed;

% parameters for the star field
STARFIELD.dimensionX = 400*coordinateMuilty;  % cm
STARFIELD.dimensionY = 400*coordinateMuilty;  % cm
STARFIELD.dimensionZ = 700*coordinateMuilty;  % cm
STARFIELD.starSize = 0.1;    % degree
STARFIELD.density = 1000/(100*coordinateMuilty)^3;    % convert num/m^3 to num/cm^3
STARFIELD.probability = TRIALINFO.coherence;

% parameters for the camera
CAMERA.elevation = 0*coordinateMuilty; % unit cm, average height of a human
CAMERA.distance = SCREEN.distance; % unit cm, distance from participant to the screen
CAMERA.sightDegreeVer = atand(SCREEN.heightCM * 0.5 / CAMERA.distance)*2; % degree of view field on vertical
CAMERA.sightDegreeHor = atand(SCREEN.widthCM * 0.5 / CAMERA.distance)*2; % degree of view field on horizon

global GL;
if testMode
    Screen('Preference', 'SkipSyncTests', 1); % for debug/test
else
    Screen('Preference', 'SkipSyncTests', 0); % for recording
end

AssertOpenGL;
InitializeMatlabOpenGL;

SCREEN.screenId = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');

% Define background color:
whiteBackground = WhiteIndex(SCREEN.screenId);
blackBackground = BlackIndex(SCREEN.screenId);

% Open a double-buffered full-screen window on the main displays screen.
[win , winRect] = PsychImaging('OpenWindow', SCREEN.screenId, whiteBackground);
SCREEN.widthPix = winRect(3);
SCREEN.heightPix = winRect(4);
SCREEN.center = [SCREEN.widthPix/2, SCREEN.heightPix/2];

TRIALINFO.fixationSizeP = degree2pix(TRIALINFO.fixationSizeD/2);
TRIALINFO.fixationPosition{1} = [SCREEN.widthPix/2-degree2pix(TRIALINFO.fixationDegree,1),SCREEN.heightPix/2];
TRIALINFO.fixationPosition{2} = [SCREEN.widthPix/2,SCREEN.heightPix/2];
TRIALINFO.fixationPosition{3} = [SCREEN.widthPix/2+degree2pix(TRIALINFO.fixationDegree,1),SCREEN.heightPix/2];

SCREEN.refreshRate = Screen('NominalFrameRate', SCREEN.screenId);

%% the configuration of the Frustum
FRUSTUM.clipNear = SCREEN.distance; % cm
FRUSTUM.clipFar = 150*coordinateMuilty; % cm
FRUSTUM.top = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.heightCM / 2.0);
FRUSTUM.bottom = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.heightCM / 2.0);
FRUSTUM.right = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 );
FRUSTUM.left = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 );


Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
glMatrixMode(GL.PROJECTION);
glLoadIdentity;
glFrustum( FRUSTUM.left,FRUSTUM.right, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);

Screen('EndOpenGL', win);

%% trial conditions and order
[TRIALINFO.trialConditions,conditionIndex] = calculateConditions();
trialIndex = repmat(conditionIndex, TRIALINFO.repetition,1);
trialNum = size(trialIndex,1);
trialOrder = randperm(trialNum);

GenerateStarField();

% calculate for the position of fixation point
[fixX,fixY] = calculateFixation();

% initial Eyelink
timePredicted = (TRIALINFO.fixationPeriod + TRIALINFO.pausePeriod + TRIALINFO.preMoveDuration+TRIALINFO.moveDuration + TRIALINFO.intertrialInterval) * TRIALINFO.repetition * length(TRIALINFO.motionType) * length(TRIALINFO.rotationDegree);
calibrationInterval = 600; % unit second, it is better to re-calibration every 10-15 minutes
automaticCalibration = timePredicted > 1.3*calibrationInterval; % make automatic calibration (every 10 min in default) if the block takes more than 15 min.

if ~testMode
    tempName = 'TEMP1'; % need temp name because Eyelink only know hows to save names with 8 chars or less. Will change name using matlab's moveFile later.
    dummymode=0;
    
    el=EyelinkInitDefaults(win);
%     el.backgroundcolour = backgroundColor;
%     el.foregroundcolour = BlackIndex(el.window);
%     el.msgfontcolour    = BlackIndex(el.window);
%     el.imgtitlecolour   = BlackIndex(el.window);
    
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        Eyelink('ShutDown');
        Screen('CloseAll');
        return
    end
    
    trialI = Eyelink('Openfile', tempName);
    if trialI~=0
        fprintf('Cannot create EDF file ''%s'' ', fileName);
        cleanup;
        Eyelink('ShutDown');
        Screen('CloseAll');
        return
    end
    
    %   SET UP TRACKER CONFIGURATION
    Eyelink('command', 'calibration_type = HV9');
    %	set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');
    Eyelink('command', 'online_dcorr_refposn = %1d, %1d', SCREEN.center(1), SCREEN.center(2));
    Eyelink('command', 'online_dcorr_maxangle = %1d', 30.0);
    % you must call this function to apply the changes from above
    EyelinkUpdateDefaults(el);
    
    % Calibrate the eye tracker
    EyelinkDoTrackerSetup(el);
    
    % do a final check of calibration using driftcorrection
    EyelinkDoDriftCorrection(el);
    
    Eyelink('StartRecording');
    
    Eyelink('message', 'SYNCTIME');	 	 % zero-plot time for EDFVIEW
    eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
    if eye_used == el.BINOCULAR % if both eyes are tracked
        eye_used = el.LEFTEYE; % use left eye
    end
    errorCheck=Eyelink('checkrecording'); 		% Check recording status */
    if(errorCheck~=0)
        fprintf('Eyelink checked wrong status.\n');
        cleanup;  % cleanup function
        Eyelink('ShutDown');
        Screen('CloseAll');
    end
    
    calibrateCkeck = tic;
    WaitSecs(1); % wait a little bit, in case the key press during calibration influence the following keyboard check
end

trialStTime = zeros(trialNum,1);
blockSt = tic;
trialI = 1;
retryFlag = 0;

while trialI <= trialNum
    
    [ ~, ~, keyCode] = KbCheck;
    if keyCode(escape)
        break;
    end
    if ~testMode
        % auto-calibration, force calibration is coded on fixation check period
        if automaticCalibration
            if toc(calibrateCkeck) >= calibrationInterval
                EyelinkDoTrackerSetup(el);
                Eyelink('StartRecording');
                Eyelink('message', 'Calibrate Finished');
                errorCheck=Eyelink('checkrecording'); 		% Check recording status */
                if(errorCheck~=0)
                    fprintf('Eyelink checked wrong status.\n');
                    cleanup;  % cleanup function
                    Eyelink('ShutDown');
                    Screen('CloseAll');
                end
                calibrateCkeck = tic;
                WaitSecs(1); % wait a little bit, in case the key press during calibration influence the following keyboard check
            end
        end
    end
    
    motionTypeI = trialIndex(trialOrder(trialI),1);
    
    % White during the inter-trial intervals
    Screen('FillRect', win ,whiteBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
    Screen('Flip', win,0,0);
    trialInterval = tic;
    
    % calculate for pre-movement
    [pglX,pglY,pglZ,pfX,pfY,pfZ] = calculatePreMove();
    
    % calculate for movement
    [glX,glY,glZ,fX,fY,fZ] = calculateMovement(motionTypeI,trialIndex(trialOrder(trialI),:),pglX(end),pglY(end),pglZ(end));
    
    % fixationType 1: L to R, 2: stay at center, 3: R to L
    if motionTypeI == 1
        fixationType = 2;
    elseif motionTypeI == 3
        fixationType = 2;
    else
        fixationType = TRIALINFO.trialConditions{trialIndex(trialOrder(trialI),1)}(trialIndex(trialOrder(trialI),2),4);
    end
    
    % wait for trial interval
    WaitSecs(TRIALINFO.intertrialInterval-toc(trialInterval)); % ITI
    
    Screen('FillRect', win ,blackBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
    drawFixation(TRIALINFO.fixationPosition{fixationType},TRIALINFO.fixationSizeP,win);
    Screen('Flip', win); % T(-2)
    
    if ~testMode
        % fixation check
        [escFlag,retryFlag] = fixationCheck(TRIALINFO.fixationPosition{fixationType},degree2pix(TRIALINFO.fixationWindow),TRIALINFO.fixationPeriod,escape,skipKey,cKey,el);
        Eyelink('message', ['Moving Start ' num2str(trialI)]);
        trialStTime(trialI) = toc(blockSt);
    else
        trialStTime(trialI) = toc(blockSt);
        escFlag=0;
        retryFlag=0;
    end
    
    if escFlag
        break;
    end
    
    if ~retryFlag
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        
        % for pre-movement
        for f=1:length(pglX)
            
            if(mod(f,1)==0)
                modifyStarField();
            end
            Screen('BeginOpenGL', win);
            glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
            glMatrixMode(GL.PROJECTION);
            glLoadIdentity;
            glFrustum( FRUSTUM.left,FRUSTUM.right, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;
            gluLookAt(pglX(f),pglY(f),pglZ(f),pfX(f),pfY(f),pfZ(f),0.0,1.0,0.0);
            glClearColor(0,0,0,0);
            glColor3f(1,1,0);
            Screen('EndOpenGL', win);
            
            % draw the fixation point and 3d dots
            DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
            drawFixation(TRIALINFO.fixationPosition{fixationType},TRIALINFO.fixationSizeP,win);
            
            % check for eye pursuit
            if ~testMode
                retryFlag = pursuitCheck(TRIALINFO.fixationPosition{fixationType},degree2pix(TRIALINFO.pursuitWindow));
            end
            if retryFlag
                break
            end
            Screen('Flip', win, 0, 0);
        end
        
        if ~retryFlag
            % for trial movement
            if ~testMode
                eyePIndex = nan(10,2);
                eyePT = tic;
                eyePi = 1;
                
                % extract current eye position
                evt = Eyelink( 'NewestFloatSample');
                eyeUsed = Eyelink('EyeAvailable'); % get eye that's tracked
                if eyeUsed ~= -1 % do we know which eye to use yet?
                    px =evt.gx(eyeUsed+1); % +1 as we're accessing MATLAB array
                    py = evt.gy(eyeUsed+1);
                    % frameStartTime(i) = evt.time;
                end
                eyePO = [px,py];
            else
                eyePO = [fixX{fixationType}(1),fixY{fixationType}(1)];
            end
            eyePNew = [fixX{fixationType}(1),fixY{fixationType}(1)];
            % define origin face direction for motion type 4
            if motionTypeI ==4
                faceDirection = [0;0;-1];
            end
            
            for f=1:length(glX)
                if(mod(f,1)==0)
                    modifyStarField();
                end
                Screen('BeginOpenGL', win);
                glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
                glMatrixMode(GL.PROJECTION);
                glLoadIdentity;
                glFrustum( FRUSTUM.left,FRUSTUM.right, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
                glMatrixMode(GL.MODELVIEW);
                glLoadIdentity;
                if motionTypeI ~= 4
                    gluLookAt(glX(f),glY(f),glZ(f),fX(f),fY(f),fZ(f),0.0,1.0,0.0);
                elseif motionTypeI == 4
                    if ~testMode
                        for k = 1:10
                            if toc(eyePT) < 0.01 && isnan(eyePIndex(end,1)) % 10ms
                                evt = Eyelink( 'NewestFloatSample');
                                eyeUsed = Eyelink('EyeAvailable'); % get eye that's tracked
                                if eyeUsed ~= -1 % do we know which eye to use yet?
                                    px =evt.gx(eyeUsed+1); % +1 as we're accessing MATLAB array
                                    py = evt.gy(eyeUsed+1);
                                    % frameStartTime(i) = evt.time;
                                end
                                eyePIndex(eyePi,:) = [px,py];
                                eyePi = eyePi+1;
                            else
                                % calculate mean position
                                eyePNew = nanmean(eyePIndex,1);
                                
                                eyeO2C = eyePO - SCREEN.center;
                                eyeN2C = eyePNew - SCREEN.center;
                                
                                % calculate for rotation on
                                eyeRD = (pix2degree(eyeN2C) - pix2degree(eyeO2C));
                                faceDirection = roty(eyeRD(1)) * (rotx(eyeRD(2))*faceDirection);
                                eyePO = eyePNew;
                                eyePT = tic;
                                eyePi = 1;
                                eyePIndex = nan(10,2);
                            end
                        end
                    else
                        eyePNew = [fixX{fixationType}(f),fixY{fixationType}(f)];
                        eyeO2C = eyePO - SCREEN.center;
                        eyeN2C = eyePNew - SCREEN.center;
                        
                        % calculate for rotation on
                        eyeRD = (pix2degree(eyeN2C) - pix2degree(eyeO2C));
                        faceDirection = roty(eyeRD(1)) * (rotx(eyeRD(2))*faceDirection);
                        eyePO = eyePNew;
                    end
                    gluLookAt(glX(f),glY(f),glZ(f),glX(f)+faceDirection(1),glY(f)+faceDirection(2),glZ(f)+faceDirection(3),0.0,1.0,0.0);
                end
                
                glClearColor(0,0,0,0);
                glColor3f(1,1,0);
                Screen('EndOpenGL', win);
                
                % draw the fixation point
                DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
                drawFixation([fixX{fixationType}(f),fixY{fixationType}(f)],TRIALINFO.fixationSizeP,win);
                drawFixation(eyePNew,TRIALINFO.fixationSizeP,win);
                % check for eye pursuit
                if ~testMode
                    retryFlag = pursuitCheck(TRIALINFO.fixationPosition{fixationType},degree2pix(TRIALINFO.pursuitWindow));
                end
                if retryFlag
                    break
                end
                Screen('Flip', win, 0, 0);
            end
        end
    end
    
    if retryFlag
        trialOrder = [trialOrder, trialOrder(trialI)];
        trialOrder(trialI) = [];
    else
        trialI = trialI + 1;
    end
    
end

Screen('Flip', win);

if ~testMode
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    try
        fprintf('Receiving data file ''%s''\n',fileName);
        status=Eyelink('ReceiveFile',tempName ,saveDir,1);
        if status > 0
            fprintf('ReceiveFile status %d\n ', status);
        end
        if exist(fileName, 'file')==2
            fprintf('Data file ''%s'' can be found in '' %s\n',fileName, pwd);
        end
    catch
        fprintf('Problem receiving data file ''%s''\n',fileName);
    end
    
    cd (saveDir);
    save(fullfile(saveDir,fileName));
    movefile([saveDir,'\',tempName,'.edf'],[saveDir,'\',fileName,'.edf']);
    
    % shut down the eyelink
    Eyelink('ShutDown');
end

%% save the real and the choiced heading
% save(fullfile(saveDir,fileName),'upTarget','lowerTarget','fixationPoint','trialStTime','fixationFinTime','choiceStTime','trialEndTime','trialCondition','choiceFinTime','fixDuration','SCREEN','trialDir','distractor');
Screen('CloseAll');
cd(curdir);
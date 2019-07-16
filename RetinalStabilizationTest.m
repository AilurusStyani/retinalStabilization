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

% set keyboard
KbName('UnifyKeyNames'); 
skipKey = KbName('space');
escape = KbName('f1');
leftArror = KbName('LeftArrow');
rightArror = KbName('RightArrow');
upArror = KbName('UpArrow');
cKey = KbName('c');

testMode = 1; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC

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
TRIALINFO.motionType = [1 2 3 4]; % 1: fixation; 2: normal pursuit; 3: simulated pursuit; 4:stabilized pursuit
TRIALINFO.headingDegree = 0 ; % degree, currently unused
TRIALINFO.headingSpeed = 50*coordinateMuilty; % cm/s
TRIALINFO.coherence = 100;
TRIALINFO.fixationSizeD = 0.25;  % degree

% fixation period (fixation point only)          --> 
% pause period (3D cloud appear, not moving)	--> 
% pre movement duration (heading only)          -->
% movement duration (heading with pursuit/fixation)
TRIALINFO.fixationPeriod = 0.5; % second
TRIALINFO.pausePeriod = 0.18; % second
TRIALINFO.preMoveDuration = 0.4; % second
TRIALINFO.moveDuration = 1; % second

TRIALINFO.fixationWindow = 1.75; % degree
TRIALINFO.pursuitWindow = 2; % degree

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
STARFIELD.dimensionZ = 900*coordinateMuilty;  % cm
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
    el.backgroundcolour = backgroundColor;
    el.foregroundcolour = BlackIndex(el.window);
    el.msgfontcolour    = BlackIndex(el.window);
    el.imgtitlecolour   = BlackIndex(el.window);
    
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup;  % cleanup function
        Eyelink('ShutDown');
        Screen('CloseAll');
        return
    end
    
    i = Eyelink('Openfile', tempName);
    if i~=0
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
end

trialStTime = zeros(trialNum,1);
blockSt = tic;

for i =1:trialNum
    
    [ ~, ~, keyCode] = KbCheck;
    if ~testMode
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
            end
        end
        if keyCode(cKey) % press c to calibrate
            EyelinkDoTrackerSetup(el);
            Eyelink('StartRecording');
            Eyelink('message', 'Force Calibrate Finished');
            errorCheck=Eyelink('checkrecording'); 		% Check recording status */
            if(errorCheck~=0)
                fprintf('Eyelink checked wrong status.\n');
                cleanup;  % cleanup function
                Eyelink('ShutDown');
                Screen('CloseAll');
            end
            calibrateCkeck = tic;
            WaitSecs(0.5); % wait a little bit, in case the key press during calibration influence the following keyboard check
        end
    end
    
    motionTypeI = trialIndex(trialOrder(i),1);
    
    % White during the inter-trial intervals
    Screen('FillRect', win ,whiteBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
    Screen('Flip', win,0,0);
    trialInterval = tic;
    
    % calculate for pre-movement
    [pglX,pglY,pglZ,pfX,pfY,pfZ] = calculatePreMove();
    
    % calculate for movement
    [glX,glY,glZ,fX,fY,fZ] = calculateMovement(motionTypeI,trialIndex(trialOrder(i),:),pglX(end),pglY(end),pglZ(end));
    
    % fixationType 1: L to R, 2: stay at center, 3: R to L
    if motionTypeI == 1
        fixationType = 2;
    elseif motionTypeI == 3
        fixationType = 2;
    else
        fixationType = TRIALINFO.trialConditions{trialIndex(trialOrder(i),1)}(trialIndex(trialOrder(i),2),4);
    end
    
    [ ~, ~, keyCode] = KbCheck;
    if keyCode(escape)
        break;
    end
    
    % wait for trial interval
    WaitSecs(TRIALINFO.intertrialInterval-toc(trialInterval)); % ITI
    
    
    Screen('FillRect', win ,blackBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
    drawFixation(TRIALINFO.fixationPosition{fixationType},TRIALINFO.fixationSizeP,win);
    Screen('Flip', win); % T(-2)
    
    if ~testMode
        [escFlag,retryFlag] = fixationCheck(TRIALINFO.fixationPosition{fixationType},degree2pix(TRIALINFO.fixationWindow),TRIALINFO.fixationPeriod,eye_used,escape,skipKey,cKey,el);
        if Eyelink( 'NewFloatSampleAvailable')>0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
        end
        Eyelink('message', ['Moving Start ' num2str(indexI)]);
        trialStTime(i) = evt.time;
    else
        trialStTime(i) = toc(blockSt);
    end
    
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
        
        Screen('Flip', win, 0, 0);
    end
    
    % for trial movement
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
        else
            % eyelink
            gluLookAt(glX(f),glY(f),glZ(f),fX(f),fY(f),fZ(f),0.0,1.0,0.0);
        end
        
        glClearColor(0,0,0,0);
        glColor3f(1,1,0);
        Screen('EndOpenGL', win);

        % draw the fixation point
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        drawFixation([fixX{fixationType}(f),fixY{fixationType}(f)],TRIALINFO.fixationSizeP,win);
        
        Screen('Flip', win, 0, 0);
    end
    
end

Screen('Flip', win);
Screen('CloseAll');









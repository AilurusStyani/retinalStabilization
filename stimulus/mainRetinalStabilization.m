clear all;
close all;
clc;

global TRIALINFO;
global CAMERA;
global STARFIELD;
global STARDATA;
global SCREEN;
global FRUSTUM;

%% path name and file name
subjectName = inputdlg({'Please input participant''s initials.'},'1',1);
fileName = ['retinalStabilization_' subjectName{1} '_' datestr(now,'yymmddHHMM')];
saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

% set keyboard
KbName('UnifyKeyNames'); 
skipKey = KbName('space');
escape = KbName('f1');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
upArror = KbName('UpArrow');
cKey = KbName('c'); % force calibration

pageUp = KbName('pageup'); % increase binocular deviation
pageDown = KbName('pagedown'); % decrease binocular deviation

testMode = 1; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC
feedback = 1; % in practice block, set 1 to provide feedback. otherwise set 0
feedbackDuration = 1; % unit s

TRIALINFO.deviation = 4; % initial binocular deviation, cm
deviationAdjust = 0.2; % how fast to adjust the deviation by key pressing, cm

%% parameters
coordinateMuilty = 1; % convert cm to coordinate system for moving distance etc.
SCREEN.distance = 60*coordinateMuilty;% cm

if testMode
%     SCREEN.widthCM = 34.5*coordinateMuilty; % cm, need to measure in your own PC
%     SCREEN.heightCM = 19.7*coordinateMuilty; % cm, need to measure in your own PC
    SCREEN.widthCM = 37.5*coordinateMuilty; % cm, need to measure in your own PC
    SCREEN.heightCM = 30*coordinateMuilty; % cm, need to measure in your own PC
else
    SCREEN.widthCM = 120*coordinateMuilty; % cm
    SCREEN.heightCM = 65*coordinateMuilty; % cm
end

TRIALINFO.repetition = 15;
TRIALINFO.motionType = [1 2 3 4]; % 1: fixation; 2: normal pursuit; 3: simulated pursuit; 4:stabilized pursuit
% TRIALINFO.headingDegree = [0] ; % degre
TRIALINFO.headingDegree = [-45 -30 -15 0 15 30 45]; % degree
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

TRIALINFO.fixationWindow = 2; % degree
TRIALINFO.pursuitWindow = 4; % degree

TRIALINFO.intertrialInterval = 1; % second

% for motion type 3
TRIALINFO.rotationDegree = [-10 10]; % ¡ã£¬the degree of the star' rotation
TRIALINFO.rotationSpeed = 10;  % ¡ã/s

% for motion type 2 and 4
TRIALINFO.fixMoveDirection = [1 3]; % only for motion type 2 and 4. 1: left to right; 2: constant at the center; 3: right to left;
TRIALINFO.fixationInitialDegree = max(TRIALINFO.rotationDegree/2); % degree ¡À to center
TRIALINFO.fixationDegree = TRIALINFO.fixationInitialDegree-1; % degree ¡À to center
TRIALINFO.fixSpeed = TRIALINFO.rotationSpeed;

% the f is the Head-Referenced value represented the distance from head to HREF plane in unit. It can be used to directly
% calculate the visual totation degree regardless the screen parameter and the distance from suject to screen.
fHREF = 15000;

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

% parameter for choice
choicePeriod = 2; % s

%% trial conditions and order
[TRIALINFO.trialConditions,conditionIndex] = calculateConditions();
trialIndex = repmat(conditionIndex, TRIALINFO.repetition,1);
trialNum = size(trialIndex,1);
trialOrder = randperm(trialNum);

disp(['This block has  ' num2str(trialNum) ' trials']);

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
calculateFrustum(coordinateMuilty);

GenerateStarField();

Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
% glEnable(GL_BLEND);
% glEnable(GL_ALPHA_BLEND_CORRECTLY);
% glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('EndOpenGL', win);

% calculate for the position of fixation point
[fixX,fixY] = calculateFixation();

% initial Eyelink
timePredicted = (TRIALINFO.fixationPeriod + TRIALINFO.pausePeriod + TRIALINFO.preMoveDuration+TRIALINFO.moveDuration + ...
    TRIALINFO.intertrialInterval + choicePeriod) * TRIALINFO.repetition * length(TRIALINFO.motionType) * length(TRIALINFO.rotationDegree) * length(TRIALINFO.headingDegree);
disp(['This block will cost  ' num2str(timePredicted/60) ' minutes']);
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
    
    testi = Eyelink('Openfile', tempName);
    if testi~=0
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
    Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
    Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,HREF,GAZERES,AREA,STATUS,INPUT,HTARGET');   
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

HideCursor(SCREEN.screenId);

trialStTime = zeros(trialNum,1);
blockSt = tic;
triali = 1;
retryFlag = 0;
choice = zeros(trialNum,1); % 0: no choice; 1: choice left; 2: choice right
choiceTime = zeros(trialNum,1);
Conditions = cell(trialNum,1);
while triali <= trialNum
    
    motionTypeI = trialIndex(trialOrder(triali),1);
    
    % fixationType 1: L to R, 2: stay at center, 3: R to L
    if motionTypeI == 1
        fixationType = 2;
    elseif motionTypeI == 3
        fixationType = 2;
    else
        fixationType = TRIALINFO.trialConditions{trialIndex(trialOrder(triali),1)}(trialIndex(trialOrder(triali),2),4);
    end
    
    % reset mouse position
    if testMode
        if motionTypeI == 4
            ShowCursor('Arrow',SCREEN.screenId);
            SetMouse(fixX{fixationType}(1),fixY{fixationType}(1))
        else
            HideCursor(SCREEN.screenId);
        end
    else
        HideCursor(SCREEN.screenId);
    end
    
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
    
    % White during the inter-trial intervals
    Screen('FillRect', win ,whiteBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
    Screen('Flip', win,0,0);
    trialInterval = tic;
    
    % calculate for pre-movement
    [pglX,pglY,pglZ,pfX,pfY,pfZ,pthetaY] = calculatePreMove(motionTypeI,trialIndex(trialOrder(triali),:));
    
    % calculate for movement
    [glX,glY,glZ,fX,fY,fZ,thetaY] = calculateMovement(motionTypeI,trialIndex(trialOrder(triali),:),pglX(end),pglY(end),pglZ(end),pfX(end),pfY(end),pfZ(end));
    
    % calculate for binocular 3D
    [pglXl,pglYl,pglZl,pfXl,pfYl,pfZl,pglXr,pglYr,pglZr,pfXr,pfYr,pfZr] = calculateCameraPosition(pglX,pglY,pglZ,pfX,pfY,pfZ);
    [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ);
    
    % wait for trial interval
    WaitSecs(TRIALINFO.intertrialInterval-toc(trialInterval)); % ITI
    
    Screen('FillRect', win ,blackBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
    Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    drawFixation(TRIALINFO.fixationPosition{fixationType},TRIALINFO.fixationSizeP,win);
    Screen('Flip', win); % T(-2)
    
    if ~testMode
        % fixation check
        escFlag = fixationCheck(TRIALINFO.fixationPosition{fixationType},degree2pix(TRIALINFO.fixationWindow),TRIALINFO.fixationPeriod,escape,skipKey,cKey,el);
        Eyelink('message', ['Moving Start ' num2str(triali)]);
        trialStTime(triali) = toc(blockSt);
    else
        trialStTime(triali) = toc(blockSt);
        escFlag = 0;
    end
    
    if escFlag
        break;
    end
    
    % for pre-movement
    for f=1:length(pglX)
        
%         if(mod(f,1)==0)
%             modifyStarField();
%         end
        
        % adjust binocular deviation
        [ ~, ~, keyCode] = KbCheck;
        if keyCode(pageUp)
            TRIALINFO.deviation = TRIALINFO.deviation + deviationAdjust;
            disp(['binocular deviation: ' num2str(TRIALINFO.deviation)]);
            [pglXl,pglYl,pglZl,pfXl,pfYl,pfZl,pglXr,pglYr,pglZr,pfXr,pfYr,pfZr] = calculateCameraPosition(pglX,pglY,pglZ,pfX,pfY,pfZ);
            [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ);
        elseif keyCode(pageDown)
            if TRIALINFO.deviation > deviationAdjust
                TRIALINFO.deviation = TRIALINFO.deviation - deviationAdjust;
                disp(['binocular deviation: ' num2str(TRIALINFO.deviation)]);
                [pglXl,pglYl,pglZl,pfXl,pfYl,pfZl,pglXr,pglYr,pglZr,pfXr,pfYr,pfZr] = calculateCameraPosition(pglX,pglY,pglZ,pfX,pfY,pfZ);
                [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ);
            end
        end
        
        calculateFrustum(coordinateMuilty,pthetaY(f));

       %% draw for left eye
        Screen('BeginOpenGL', win);
        glColorMask(GL.TRUE, GL.FALSE, GL.FALSE, GL.FALSE);
        glMatrixMode(GL.PROJECTION);
        glLoadIdentity;
        glFrustum( FRUSTUM.sinisterLeft,FRUSTUM.sinisterRight, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
        glMatrixMode(GL.MODELVIEW);
        glLoadIdentity;
        gluLookAt(pglXl(f),pglYl(f),pglZl(f),pglXl(f)+pfXl(f),pglYl(f)+pfYl(f),pglZl(f)+pfZl(f),0.0,1.0,0.0);
        glClearColor(0,0,0,0);
        glColor3f(1,1,0);
        
        % draw the fixation point and 3d dots
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        
        %% draw for right eye
        glColorMask(GL.FALSE, GL.TRUE, GL.FALSE, GL.FALSE);
        glMatrixMode(GL.PROJECTION);
        glLoadIdentity;
        glFrustum( FRUSTUM.dexterLeft,FRUSTUM.dexterRight, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
        glMatrixMode(GL.MODELVIEW);
        glLoadIdentity;
        gluLookAt(pglXr(f),pglYr(f),pglZr(f),pglXr(f)+pfXr(f),pglYr(f)+pfYr(f),pglZr(f)+pfZr(f),0.0,1.0,0.0);
        glClearColor(0,0,0,0);
        glColor3f(1,1,0);
        
        
        % draw the fixation point and 3d dots for right eye
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        Screen('EndOpenGL', win);
        
        drawFixation(TRIALINFO.fixationPosition{fixationType},TRIALINFO.fixationSizeP,win);
       %% check for eye pursuit
        if ~testMode
            retryFlag = pursuitCheck(TRIALINFO.fixationPosition{fixationType},degree2pix(TRIALINFO.pursuitWindow),win);
        end
        if retryFlag
            break
        end
        Screen('Flip', win, 0, 0);
        if f==1
            WaitSecs(TRIALINFO.pausePeriod);
        end
    end
        
    if ~retryFlag
        % for trial movement
        if ~testMode
            
            % extract current eye position
            evt = Eyelink( 'NewestFloatSample');
            eyeUsed = Eyelink('EyeAvailable'); % get eye that's tracked
            if eyeUsed ~= -1 % do we know which eye to use yet?
                hx =evt.hx(eyeUsed+1); % +1 as we're accessing MATLAB array
                hy = evt.hy(eyeUsed+1);
                % frameStartTime(i) = evt.time;
            end
            eyePO = [hx,hy];
        else
            eyePO = [fixX{fixationType}(1),fixY{fixationType}(1)];
        end
        
        % define origin face direction for motion type 4
        if motionTypeI ==4
            faceDirectionL = [fXl(1);fYl(1);fZl(1)];
            faceDirectionR = [fXr(1);fYr(1);fZr(1)];
        end
        
        for f=1:length(glX)
%             if(mod(f,1)==0)
%                 modifyStarField();
%             end
            
            % adjust binocular deviation
            [ ~, ~, keyCode] = KbCheck;
            if keyCode(pageUp)
                TRIALINFO.deviation = TRIALINFO.deviation + deviationAdjust;
                disp(['binocular deviation: ' num2str(TRIALINFO.deviation)]);
                [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ);
            elseif keyCode(pageDown)
                if TRIALINFO.deviation > deviationAdjust
                    TRIALINFO.deviation = TRIALINFO.deviation - deviationAdjust;
                    disp(['binocular deviation: ' num2str(TRIALINFO.deviation)]);
                    [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ);
                end
            end
            
            calculateFrustum(coordinateMuilty,thetaY(f));

           %% draw for left eye
            Screen('BeginOpenGL', win);
            glColorMask(GL.TRUE, GL.FALSE, GL.FALSE, GL.FALSE);
            glMatrixMode(GL.PROJECTION);
            glLoadIdentity;
            glFrustum( FRUSTUM.sinisterLeft,FRUSTUM.sinisterRight, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;
            
            if motionTypeI ~= 4
                gluLookAt(pXl(f),pYl(f),pZl(f),pXl(f)+fXl(f),pYl(f)+fYl(f),pZl(f)+fZl(f),0.0,1.0,0.0);
            elseif motionTypeI == 4
                if ~testMode
                    eyePT = tic;
                    eyePIndex = nan(10,2);
                    for k = 1:11
                        if toc(eyePT) < 0.010 && isnan(eyePIndex(end,1)) && (k<11) % 10ms
                            if Eyelink( 'NewFloatSampleAvailable')>0
                                % get the sample in the form of an event structure                                
                                evt = Eyelink( 'NewestFloatSample');
                                eyeUsed = Eyelink('EyeAvailable'); % get eye that's tracked
                                if eyeUsed ~= -1 % do we know which eye to use yet?
                                    hx =evt.hx(eyeUsed+1); % +1 as we're accessing MATLAB array
                                    hy = evt.hy(eyeUsed+1);
                                    % frameStartTime(i) = evt.time;
                                end
                                WaitSecs(0.001);
                            end
                            eyePIndex(k,:) = [hx,hy];
                            disp(['Sample position ' num2str([hx hy])])
                        else
                            % calculate mean position
                            eyePNew = nanmean(eyePIndex,1);
                            if ~isnan(eyePNew)
                                
                                % calculate for rotation degree
                                eyeRD = [atand(eyePNew(1)/fHREF)-atand(eyePO(1)/fHREF),atand(eyePNew(2)/fHREF)-atand(eyePO(2)/fHREF)];
                                
                                faceDirectionL = roty(eyeRD(1)) * (rotx(eyeRD(2))*faceDirectionL);
                                faceDirectionR = roty(eyeRD(1)) * (rotx(eyeRD(2))*faceDirectionR);
                                eyePO = eyePNew;
                            end
                            break
                        end
                    end
                else
                    [eyePNew(1),eyePNew(2),~] = GetMouse;
                    eyeO2C = eyePO - SCREEN.center;
                    eyeN2C = eyePNew - SCREEN.center;
                    
                    % calculate for rotation on
                    eyeRD = (pix2degree(eyeN2C) - pix2degree(eyeO2C));
                    faceDirectionL = roty(eyeRD(1)) * (rotx(eyeRD(2))*faceDirectionL);
                    faceDirectionR = roty(eyeRD(1)) * (rotx(eyeRD(2))*faceDirectionR);
                    eyePO = eyePNew;
                end
                gluLookAt(pXl(f),pYl(f),pZl(f),pXl(f)+faceDirectionL(1),pYl(f)+faceDirectionL(2),pZl(f)+faceDirectionL(3),0.0,1.0,0.0);
            end
            
            glClearColor(0,0,0,0);
            glColor3f(1,1,0);
            
            % draw the fixation point
            DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
            
          %% draw for right eye
            glColorMask(GL.FALSE, GL.TRUE, GL.FALSE, GL.FALSE);
            glMatrixMode(GL.PROJECTION);
            glLoadIdentity;
            glFrustum( FRUSTUM.dexterLeft,FRUSTUM.dexterRight, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
            glMatrixMode(GL.MODELVIEW);
            glLoadIdentity;
            
            if motionTypeI ~= 4
                gluLookAt(pXr(f),pYr(f),pZr(f),pXr(f)+fXr(f),pYr(f)+fYr(f),pZr(f)+fZr(f),0.0,1.0,0.0);
            elseif motionTypeI == 4
                gluLookAt(pXr(f),pYr(f),pZr(f),pXr(f)+faceDirectionR(1),pYr(f)+faceDirectionR(2),pZr(f)+faceDirectionR(3),0.0,1.0,0.0);
            end
            
            glClearColor(0,0,0,0);
            glColor3f(1,1,0);
            DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
            Screen('EndOpenGL', win);
            
            % draw the fixation point
            drawFixation([fixX{fixationType}(f),fixY{fixationType}(f)],TRIALINFO.fixationSizeP,win);
            
            % check for eye pursuit
            if ~testMode
                retryFlag = pursuitCheck([fixX{fixationType}(f),fixY{fixationType}(f)],degree2pix(TRIALINFO.pursuitWindow),win);
            else
                % simulate eye movement
%                 drawFixation(eyePNew,TRIALINFO.fixationSizeP,win);
            end
            
            if retryFlag
                break
            end
            Screen('Flip', win, 0, 0);
        end
    end
    
    if ~retryFlag
        % starting choice
        degree = TRIALINFO.trialConditions{trialIndex(trialOrder(triali),1)}(trialIndex(trialOrder(triali),2),2);
        correctAnswer = (degree >= 0)+1;
        if ~testMode
            Eyelink('message', ['Start choice ' num2str(triali)]);
        end
        choice(triali) = 0;
        startChoice = tic;
        [~, ~, ~] = DrawFormattedText(win, 'What''s your heading direction?','center',SCREEN.center(2)/2,[200 200 200]);
        Screen('TextBackgroundColor',win, [0 0 0 0]);
        Screen('DrawingFinished',win);
        Screen('Flip',win,0,0);
        while toc(startChoice) <= choicePeriod
            [ ~, ~, keyCode ] = KbCheck;
            if keyCode(leftKey)
                choice(triali) = 1;
                choiceTime(triali) = toc(startChoice);
            elseif keyCode(rightKey)
                choice(triali) = 2;
                choiceTime(triali) = toc(startChoice);
            end
            if choice(triali)
                break
            end
        end
        if feedback
            if choice(triali) == correctAnswer
                sound(sin(2*pi*25*(1:3000)/200)); % correct cue
                [~, ~, ~] = DrawFormattedText(win, 'You are right!','center',SCREEN.center(2)/2,[20 200 20]);
                if ~testMode
                Eyelink('message', ['Decision made ' num2str(triali)]);
                end
            elseif choice(triali)
                sound(sin(2*pi*25*(1:3000)/600)); % wrong cue
                [~, ~, ~] = DrawFormattedText(win, 'Please try again.','center',SCREEN.center(2)/2,[200 20 20]);
                if ~testMode
                Eyelink('message', ['Decision made ' num2str(triali)]);
                end
            else
                sound(sin(2*pi*25*(1:3000)/600)); % missing cue
                [~, ~, ~] = DrawFormattedText(win, 'Oops, you missing this trial.','center',SCREEN.center(2)/2,[200 20 20]);
                if ~testMode
                    Eyelink('message', ['Missing ' num2str(triali)]);
                end
            end
            Screen('TextBackgroundColor',win, [0 0 0 0]);
            Screen('DrawingFinished',win);
            Screen('Flip',win,0,0);
            WaitSecs(feedbackDuration);
        else
            if choice(triali)
                sound(sin(2*pi*25*(1:3000)/200)); % response cue
                if ~testMode
                    Eyelink('message', ['Decision made ' num2str(triali)]);
                end
            else
                sound(sin(2*pi*25*(1:3000)/600)); % missing cue
                if ~testMode
                    Eyelink('message', ['Missing ' num2str(triali)]);
                end
            end
        end
    end
    
    if retryFlag || ~choice(triali)
        if ~testMode
            Eyelink('message', ['Trial repeat ' num2str(triali)]);
        end
        trialOrder = [trialOrder, trialOrder(triali)];
        trialOrder(triali) = [];
    else
        Conditions{triali} = TRIALINFO.trialConditions{trialIndex(trialOrder(triali),1)}(trialIndex(trialOrder(triali),2),:);
        if ~testMode
            Eyelink('message', ['Trial complete ' num2str(triali)]);
        end
        triali = triali + 1;
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
save(fullfile(saveDir,fileName),'TRIALINFO','Conditions','SCREEN','CAMERA','choice','choiceTime');
Screen('CloseAll');
cd(curdir);
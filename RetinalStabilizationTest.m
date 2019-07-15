clear all;
close all;
clc;

global TRIALINFO;
global CAMERA;
global STARFIELD;
global STARDATA;
global SCREEN

% set keyboard
KbName('UnifyKeyNames'); 
skipKey = KbName('space');
escape = KbName('f1');
leftArror = KbName('LeftArrow');
rightArror = KbName('RightArrow');
upArror = KbName('UpArrow');
cKey = KbName('c');

testMode = 1; % in test mode, the codes related to Eyelink will be skipped so that you can debug in your own PC
calibrationInterval = 600; % unit second, it is better to re-calibration every 10-15 minutes

%% parameters
coordinateMuilty = 5; % convert cm to coordinate system for moving distance
SCREEN.distance = 60*coordinateMuilty;% cm

if testMode
    SCREEN.widthCM= 37.5*coordinateMuilty; % cm
    SCREEN.heightCM = 30*coordinateMuilty; % cm
else
    SCREEN.widthCM= 120*coordinateMuilty; % cm
    SCREEN.heightCM = 65*coordinateMuilty; % cm
end

TRIALINFO.repetition =3;
TRIALINFO.motionType = [1 2 3 4]; % 1: fixation; 2: normal pursuit; 3: simulated pursuit; 4:stabilized pursuit
TRIALINFO.preMoveDuration = 0.2; % second
TRIALINFO.moveDuration = 1; % second
TRIALINFO.headingDegree = 0 ; % degree
TRIALINFO.headingSpeed = 50*coordinateMuilty; % cm/s
TRIALINFO.coherence = 100;
TRIALINFO.fixationSizeD = 0.25;  % degree
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
STARFIELD.dimensionZ = 400*coordinateMuilty;  % cm
STARFIELD.starSize = 0.1;    % degree
STARFIELD.density = 1000/(100*coordinateMuilty)^3;    % convert num/m^3 to num/cm^3
STARFIELD.probability = TRIALINFO.coherence;

% parameters for the camera
CAMERA.elevation = 0*coordinateMuilty; % unit cm, average height of a human
CAMERA.distance = SCREEN.distance; % unit cm, distance from participant to the screen
CAMERA.sightDegreeVer = atand(SCREEN.heightCM * 0.5 / CAMERA.distance)*2; % degree of view field on vertical
CAMERA.sightDegreeHor = atand(SCREEN.widthCM * 0.5 / CAMERA.distance)*2; % degree of view field on horizon

global GL;
AssertOpenGL;
InitializeMatlabOpenGL;

SCREEN.screenId = max(Screen('Screens'))-1;
PsychImaging('PrepareConfiguration');

if testMode
    Screen('Preference', 'SkipSyncTests', 1); % for debug/test
else
    Screen('Preference', 'SkipSyncTests', 0); % for recording
end

% Open a double-buffered full-screen window on the main displays screen.
[win, winRect] = PsychImaging('OpenWindow', SCREEN.screenId, 0, [], [], [], 0, 0);
SCREEN.widthPix = winRect(3);
SCREEN.heightPix = winRect(4);
SCREEN.pixInDegree = SCREEN.widthPix/CAMERA.sightDegreeVer;
SCREEN.cmInDegree = SCREEN.widthCM/CAMERA.sightDegreeVer;

TRIALINFO.fixationSize = degree2pix(TRIALINFO.fixationSizeD); % pixel
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

for i =1:trialNum
    motionTypeI = trialIndex(trialOrder(i),1);
    
    % White during the inter-trial intervals
    Screen('FillRect', win ,[0 0 0],[0 0 SCREEN.widthPix SCREEN.heightPix]);
    Screen('Flip', win,0,0);
    trialInterval = tic;
    
    % calculate for pre-movement
    [pglX,pglY,pglZ,pfX,pfY,pfZ] = calculatePreMove();
    
    % calculate for movement
    [glX,glY,glZ,fX,fY,fZ] = calculateMovement(motionTypeI);
    
    % fixationType 1: L to R, 2: stay at center, 3: R to L
    if motionTypeI == 1
        fixationType = 2;
    elseif motionTypeI == 3
        fixationType = 2;
    else
        fixationType = TRIALINFO.trialConditions{trialIndex(trialOrder(i),1)}(trialIndex(trialOrder(i),2),4);
    end
    
    drawFixation(TRIALINFO.fixationPosition{fixationType});
    
    [ ~, ~, keyCode] = KbCheck;
    if keyCode(escape)
        break;
    end
    
    % wait for trial interval
    WaitSecs(TRIALINFO.intertrialInterval-toc(trialInterval)); % ITI
    
    
    Screen('Flip', win); % T(-2)
    
    %eyelink check for fixation
    
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
        glMatrixMode(GL.MODELVIEW);
        gluLookAt(pglX(f),pglY(f),pglZ(f),pfX(f),pfY(f),pfZ(f),0.0,1.0,0.0);
        glClearColor(0,0,0,0);
        glColor3f(1,1,0);
        Screen('EndOpenGL', win);
        
        % draw the fixation point
        drawFixation(TRIALINFO.fixationPosition{fixationType});
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        Screen('Flip', win, 0, 0);
    end

end

Screen('Flip', win);
Screen('CloseAll');









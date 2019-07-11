clear all;
close all;
clc;

global TRIALINFO;
global CAMERA;
global STARFIELD;
global STARDATA;
global SCREEN

% Restrict KbCheck to checking of ESCAPE key:
KbName('UnifyKeynames');
escape = KbName('ESCAPE');

%% trial information
coordinateMuilty = 5; % convert cm to coordinate system for moving distance
SCREEN.distance = 60;% cm
SCREEN.widthCM= 120; % cm
SCREEN.heightCM = 65; % cm

TRIALINFO.repeatNum =3;
TRIALINFO.motionType = [1 2 3 4];
TRIALINFO.preMoveDuration = 0.2; % s
TRIALINFO.moveDuration = 1; % second
TRIALINFO.headingDegree = 0 ; % degree
TRIALINFO.headingSpeed = 50*coordinateMuilty; % cm/s

TRIALINFO.rotationDegree = 30; % the degree of the star' rotation
TRIALINFO.rotationSpeed = 10;  % бу/s
TRIALINFO.coherence = 100;

TRIALINFO.fixMoveDirection = [-1 1];
TRIALINFO.fixzationSizeD = 0.25;  % degree
TRIALINFO.fixSpeed = TRIALINFO.rotationSpeed;

%% the information of the stars
STARFIELD.demensionX = 400*coordinateMuilty;  % cm
STARFIELD.demensionY = 400*coordinateMuilty;  % cm
STARFIELD.demensionZ = 40*coordinateMuilty;   % cm
STARFIELD.StarSize = 0.2;    % cm
STARFIELD.Density = 0.01;    % num/cm^3
STARFIELD.Probability = TRIALINFO.coherence;

%% parameter for camera
CAMERA.elevation = 0; % unit cm, average height of a human
CAMERA.distance = SCREEN.distance; % unit cm, distance from participant to the screen
CAMERA.sightDegreeVer = atand(SCREEN.heightCM * 0.5 / CAMERA.distance)*2; % degree of view field on vertical
CAMERA.sightDegreeHor = atand(SCREEN.widthCM * 0.5 / CAMERA.distance)*2; % degree of view field on horizon

global GL;
AssertOpenGL;
InitializeMatlabOpenGL;

SCREEN.screenId = max(Screen('Screens'))-1;
PsychImaging('PrepareConfiguration');

Screen('Preference', 'SkipSyncTests', 1); % for debug/test
% Screen('Preference', 'SkipSyncTests', 0); % for experiment

% Open a double-buffered full-screen window on the main displays screen.
[win, winRect] = PsychImaging('OpenWindow', SCREEN.screenId, 0, [], [], [], 0, 0);
SCREEN.widthPix = winRect(3);
SCREEN.heightPix = winRect(4);
SCREEN.pixInDegree = SCREEN.screenWidthPix/CAMERA.sightDegreeVer;
SCREEN.cmInDegree = SCREEN.screenWidthCM/CAMERA.sightDegreeVer;

TRIALINFO.fixationSize = degree2pix(TRIALINFO.fixationSizeD); % pixel
TRIALINFO.fix_move_distance = degree2pix(TRIALINFO.rotationSpeed * TRIALINFO.moveDuration);  % pixel

SCREEN.refreshRate = Screen('NominalFrameRate', screenId);

%% the configuration of the Frustum
FRUSTUM.clipNear = SCREEN.distance; % cm
FRUSTUM.clipFar = 150; % cm
FRUSTUM.top = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.heightCM / 2.0);
FRUSTUM.bottom = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.heightCM / 2.0);
FRUSTUM.right = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 );
FRUSTUM.left = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 );


Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
glFrustum( FRUSTUM.left,FRUSTUM.right, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
Screen('EndOpenGL', win);

%% trial num and order
conditions = calculateConditions();
if(TRIALINFO.motionType<3) %% fixation point is stable and the starts are moving directly ||  %% fixation point is moving directly and the starts are moving directly
    TRIALINFO.trial_num = length(TRIALINFO.headingDegree)*TRIALINFO.repeatNum;
    heading_init = repmat(TRIALINFO.headingDegree,1,TRIALINFO.repeatNum); %%% heading for star's movement
    order = randperm(TRIALINFO.trial_num);
    TRIALINFO.headings_order = heading_init(order);%% real heading order
else%% fixation point is moving directly and the starts are rotating || %% fixation point is stable and the starts are rotating
    
    TRIALINFO.trial_num = length(TRIALINFO.headingDegree)*TRIALINFO.repeatNum;
    heading_init = repmat(TRIALINFO.headingDegree,1,TRIALINFO.repeatNum); %%% heading for star's movement
    order = randperm(TRIALINFO.trial_num);
    TRIALINFO.headings_order = heading_init(order);%% real heading order
    
    rotation_init = repmat(TRIALINFO.rotationDegree,1,TRIALINFO.repeatNum); %%% rotation degree for star's movement
    TRIALINFO.rotationDegree_order = rotation_init(order);%% real heading order
end

fixation_position = [TRIALINFO.screenWidthPIX/2-TRIALINFO.fixzationSize TRIALINFO.screenHeightPIX/2-TRIALINFO.fixzationSize TRIALINFO.screenWidthPIX/2+TRIALINFO.fixzationSize TRIALINFO.screenHeightPIX/2+TRIALINFO.fixzationSize];

GenerateStarFiled(STARFIELD.demensionX,STARFIELD.demensionY,STARFIELD.demensionZ,STARFIELD.StarSize,STARFIELD.Density);
fixX = [];

for i =1:TRIALINFO.trial_num
    
    [ keyIsDown, seconds, keyCode ] = KbCheck;
    if keyCode(escape)
        break;
    end
    
     Screen('FillRect', win ,[0 0 0],[0 0 TRIALINFO.screenWidthPIX TRIALINFO.screenHeightPIX]);
  
    trialCondition = zeros(1,5);
    if(TRIALINFO.motionType<3)
        trialCondition(1) = TRIALINFO.headings_order(i);
    else
        trialCondition(1) = TRIALINFO.headings_order(i);
        trialCondition(2) = TRIALINFO.star_move_distance;
        trialCondition(3) = TRIALINFO.rotationDegree;
        
        trialCondition(4) = TRIALINFO.fix_move_distance;
        trialCondition(5) = TRIALINFO.moveDuration;
    end
    if( TRIALINFO.motionType == 1 ) %% fixation point is stable and the starts are moving directly
        [Lateral,Surge,Heave]= CalculateConstentVelocityMovement(0,trialCondition(1),TRIALINFO.star_move_distance,TRIALINFO.moveDuration,SCREEN.refreshRate);
    elseif(TRIALINFO.motionType == 2)%% fixation point is moving directly and the starts are moving directly
        [Lateral,Surge,Heave]= CalculateConstentVelocityMovement(0,trialCondition(1),TRIALINFO.star_move_distance,TRIALINFO.moveDuration,SCREEN.refreshRate);
        if(trialCondition(1)<90)
            towards = 1;
        else
            towards = -1;
        end
        [fixX,fixY,fixZ] = calculateLinearMovementForFix(TRIALINFO.fix_move_distance,TRIALINFO.moveDuration,towards);
    elseif(TRIALINFO.motionType == 3)%% fixation point is stable and the starts are rotating
% % %         [Lateral,Heave,Surge,fX,fY,fZ] = cameraRotation(trialCondition);
        [Lateral,Surge,Heave]= CalculateConstentVelocityMovement(0,trialCondition(1),TRIALINFO.star_move_distance,TRIALINFO.moveDuration,SCREEN.refreshRate);
       [rotationAngle,rotVector]= CalculateRotationMovement(90,0,trialCondition(3),trialCondition(5));
    else %% fixation point is moving directly and the starts are rotating
        [Lateral,Surge,Heave,fX,fY,fZ] = cameraRotation(trialCondition);
        if(trialCondition(1)<90)
            towards = 1;
        else
            towards = -1;
        end
        [fixX,fixY,fixZ] = calculateLinearMovementForFix(TRIALINFO.fix_move_distance,TRIALINFO.moveDuration,towards);
    end
   
    Screen('FillOval', win, [255 0 0 255], fixation_position);
    Screen('Flip', win);
    WaitSecs(0.5);
    Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    fixation_position_frame = fixation_position;
    for frames=1:length(Heave)
        Heavei = Heave(frames);
        Surgei = Surge(frames);
        Laterali = Lateral(frames);
       
        if(TRIALINFO.motionType >2)
            
            rotationAnglei =rotationAngle(frames);
% % %             fXi = fX(frames);
% % %             fYi = fY(frames);
% % %             fZi = fZ(frames);
        end
        
        
        if(~isempty(fixX))
            fixXi = fixX(frames);
        end
        if(mod(frames,1)==0)
            modifyStarField();
        end
        Screen('BeginOpenGL', win);
        glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
        glMatrixMode(GL.PROJECTION);
        glLoadIdentity;
               
        glMatrixMode(GL.MODELVIEW);
        glLoadIdentity;
        gluLookAt(Laterali, 0.0-Heavei, CAMERA.distance-Surgei,Laterali,......
                0.0-Heavei, CAMERA.distance-Surgei-1.0,0.0,1.0,0.0);
        if(TRIALINFO.motionType >2)
            rotCenterX = Laterali ;
            rotCenterY = -Heavei;
            rotCenterZ = CAMERA.distance-Surgei;
            glTranslatef(rotCenterX,rotCenterY,rotCenterZ);
            glRotated(rotationAnglei, rotVector(1), rotVector(2), rotVector(3));
            glTranslatef(-rotCenterX,-rotCenterY,-rotCenterZ);
        end
% % % %         if(TRIALINFO.motionType < 3 )
% % % %               gluLookAt(Laterali, 0.0-Heavei, CAMERA.distance-Surgei,Laterali,......
% % % %                 0.0-Heavei, CAMERA.distance-Surgei-1.0,0.0,1.0,0.0);
% % % %         else 
% % % %             gluLookAt(Laterali, Heavei, CAMERA.distance-Surgei,fXi,fYi,CAMERA.distance- fZi,0.0,1.0,0.0);
% % % % % %             glRotated(rotationAnglei, rotVector(1), rotVector(2), rotVector(3));
% % % %         end
        
        glClearColor(0,0,0,0);
        %%% draw the fixation point
        glColor3f(1,1,0);
        Screen('EndOpenGL', win);
        if( TRIALINFO.motionType == 2||TRIALINFO.motionType == 4 )
            fixation_position_frame = fixation_position_frame + [fixXi 0 fixXi 0];
        end
        Screen('FillOval', win, [255 0 0 255], fixation_position_frame);
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        Screen('Flip', win, 0, 0);
    end
    %% White during the inter-trial intervals  
    Screen('FillRect', win ,[255 255 255],[0 0 TRIALINFO.screenWidthPIX TRIALINFO.screenHeightPIX]);

    Screen('Flip', win,0,0);
    WaitSecs(1);
    [ keyIsDown, seconds, keyCode ] = KbCheck;
    if keyCode(escape)
        break;
    end
end

Screen('Flip', win);
Screen('CloseAll');









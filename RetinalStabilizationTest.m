clear all;
close all;
clc;

global TRIALINFO;
global CAMERA;
global STARFIELD;
global STARDATA;

% Restrict KbCheck to checking of ESCAPE key:
KbName('UnifyKeynames');
escape = KbName('ESCAPE');

%% trial information
TRIALINFO.camera2screenDist = 60;%cm
TRIALINFO.screenWidthCM= 120; %cm
TRIALINFO.screenHeightCM = 65; %cm

TRIALINFO.repeatNum =3;
TRIALINFO.motionType = 3 ;
TRIALINFO.move_time = 1; %second
TRIALINFO.star_move_heading = 0 ; %degree
TRIALINFO.star_rotation_degree = 30;%% the degree of the star' rotation
% TRIALINFO.fixMoveDirection = [-1 1];
TRIALINFO.star_coherence = 100;
TRIALINFO.fixzationSpeed =10;   %% 10��/s
TRIALINFO.star_move_speed = 50; % cm/s
TRIALINFO.fixzationSize = 0.25;  % degree



%% the information of the stars
STARFIELD.demensionX = 400;  % cm
STARFIELD.demensionY = 400;  % cm
STARFIELD.demensionZ = 40;   % cm
STARFIELD.StarSize = 0.2;    % cm
STARFIELD.Density = 0.01;    % num/cm^3
STARFIELD.Probability = TRIALINFO.star_coherence;

%% parameter for camera
CAMERA.elevation = 0; % unit cm, average height of a human
CAMERA.distance = TRIALINFO.camera2screenDist; % unit cm, distance from participant to the screen
CAMERA.sightDegreeVer = atand(TRIALINFO.screenHeightCM * 0.5 / CAMERA.distance)*2; % degree of view field on vertical
CAMERA.sightDegreeHor = atand(TRIALINFO.screenWidthCM * 0.5 / CAMERA.distance)*2; % degree of view field on horizon

global GL;
AssertOpenGL;
InitializeMatlabOpenGL;

TRIALINFO.screenId = max(Screen('Screens'))-1;
PsychImaging('PrepareConfiguration');
Screen('Preference', 'SkipSyncTests', 1);
% Open a double-buffered full-screen window on the main displays screen.
[win, winRect] = PsychImaging('OpenWindow', TRIALINFO.screenId, 0, [], [], [], 0, 0);
TRIALINFO.screenWidthPIX = winRect(3);
TRIALINFO.screenHeightPIX = winRect(4);
TRIALINFO.pixInDegree = TRIALINFO.screenWidthPIX/CAMERA.sightDegreeVer;
TRIALINFO.cmInDegree = TRIALINFO.screenWidthCM/CAMERA.sightDegreeVer;

TRIALINFO.star_move_distance = TRIALINFO.star_move_speed * TRIALINFO.move_time ;  % cm
TRIALINFO.fix_move_distance = TRIALINFO.fixzationSpeed * TRIALINFO.move_time * TRIALINFO.cmInDegree;  % cm
TRIALINFO.fixzationSize= ceil(TRIALINFO.fixzationSize*TRIALINFO.pixInDegree);

TRIALINFO.refreshRate = round(Screen('FrameRate', win ));

%% the configuration of the Frustum
FRUSTUM.clipNear = 50;
FRUSTUM.clipFar = 150;
FRUSTUM.top = (FRUSTUM.clipNear / TRIALINFO.camera2screenDist) * (TRIALINFO.screenHeightCM / 2.0);
FRUSTUM.bottom = (FRUSTUM.clipNear / TRIALINFO.camera2screenDist) * (-TRIALINFO.screenHeightCM / 2.0);
FRUSTUM.right = (FRUSTUM.clipNear / TRIALINFO.camera2screenDist) * (TRIALINFO.screenWidthCM / 2.0 );
FRUSTUM.left = (FRUSTUM.clipNear / TRIALINFO.camera2screenDist) * (-TRIALINFO.screenWidthCM / 2.0 );


Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
Screen('EndOpenGL', win);

%% trial num and order
if(TRIALINFO.motionType<3) %% fixation point is stable and the starts are moving directly ||  %% fixation point is moving directly and the starts are moving directly
    TRIALINFO.trial_num = length(TRIALINFO.star_move_heading)*TRIALINFO.repeatNum;
    heading_init = repmat(TRIALINFO.star_move_heading,1,TRIALINFO.repeatNum); %%% heading for star's movement
    order = randperm(TRIALINFO.trial_num);
    TRIALINFO.headings_order = heading_init(order);%% real heading order
else%% fixation point is moving directly and the starts are rotating || %% fixation point is stable and the starts are rotating
    
    TRIALINFO.trial_num = length(TRIALINFO.star_move_heading)*TRIALINFO.repeatNum;
    heading_init = repmat(TRIALINFO.star_move_heading,1,TRIALINFO.repeatNum); %%% heading for star's movement
    order = randperm(TRIALINFO.trial_num);
    TRIALINFO.headings_order = heading_init(order);%% real heading order
    
    rotation_init = repmat(TRIALINFO.star_rotation_degree,1,TRIALINFO.repeatNum); %%% rotation degree for star's movement
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
        trialCondition(3) = TRIALINFO.star_rotation_degree;
        
        trialCondition(4) = TRIALINFO.fix_move_distance;
        trialCondition(5) = TRIALINFO.move_time;
    end
    if( TRIALINFO.motionType == 1 ) %% fixation point is stable and the starts are moving directly
        [Lateral,Surge,Heave]= CalculateConstentVelocityMovement(0,trialCondition(1),TRIALINFO.star_move_distance,TRIALINFO.move_time);
    elseif(TRIALINFO.motionType == 2)%% fixation point is moving directly and the starts are moving directly
        [Lateral,Surge,Heave]= CalculateConstentVelocityMovement(0,trialCondition(1),TRIALINFO.star_move_distance,TRIALINFO.move_time);
        if(trialCondition(1)<90)
            towards = 1;
        else
            towards = -1;
        end
        [fixX,fixY,fixZ] = calculateLinearMovementForFix(TRIALINFO.fix_move_distance,TRIALINFO.move_time,towards);
    elseif(TRIALINFO.motionType == 3)%% fixation point is stable and the starts are rotating
% % %         [Lateral,Heave,Surge,fX,fY,fZ] = cameraRotation(trialCondition);
        [Lateral,Surge,Heave]= CalculateConstentVelocityMovement(0,trialCondition(1),TRIALINFO.star_move_distance,TRIALINFO.move_time);
       [rotationAngle,rotVector]= CalculateRotationMovement(90,0,trialCondition(3),trialCondition(5));
    else %% fixation point is moving directly and the starts are rotating
        [Lateral,Surge,Heave,fX,fY,fZ] = cameraRotation(trialCondition);
        if(trialCondition(1)<90)
            towards = 1;
        else
            towards = -1;
        end
        [fixX,fixY,fixZ] = calculateLinearMovementForFix(TRIALINFO.fix_move_distance,TRIALINFO.move_time,towards);
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
        
        glFrustum( FRUSTUM.left,FRUSTUM.right, FRUSTUM.bottom, FRUSTUM.top, FRUSTUM.clipNear, FRUSTUM.clipFar);
        
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









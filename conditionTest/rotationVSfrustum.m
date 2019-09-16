% test 3D in OPENGL
global SCREEN
global STARDATA
global STARFIELD
global MOVE

%% path name and file name
subjectName = inputdlg({'Please input participant''s initials.'},'1',1);
fileName = ['retinalStabilization_' subjectName{1} '_' datestr(now,'yymmddHHMM')];
saveDir = fullfile(pwd,'data');
mkdir(saveDir);
curdir = pwd;

testMode = 1; % 0 for record, 1 for debug
feedback = 1;
feedbackDuration = 1; % second

STARFIELD.starSize = 0.5;

repetition = 10;

deviation = 1.6; % initial binocular deviation, cm
deviationAdjust = 0.2; % how fast to adjust the deviation by key pressing, cm

%% parameters
SCREEN.distance = 60;% cm
% SCREEN.widthCM = 37.5; % cm, need to measure in your own PC
% SCREEN.heightCM = 30; % cm, need to measure in your own PC
SCREEN.widthCM = 53; % cm, need to measure in your own PC
SCREEN.heightCM = 30; % cm, need to measure in your own PC

MOVE.degree = [-10 -5 -1 1 5 10]; % бу
MOVE.velocity = 50; % cm/s
MOVE.duration = 1; % second
MOVE.rotationDegree = [-10 10]; % бу
MOVE.initialDegree = max(MOVE.rotationDegree)/2; % бу

choicePeriod = 2; % s

STARFIELD.dimensionX = 200;
STARFIELD.dimensionY = 200;
STARFIELD.dimensionZ = 200;
STARFIELD.density = 0.01;

conditions = calculateCondition();
conditions = repmat(conditions,repetition,1);
trialNum = size(conditions,1);
tempIndex = [(1:trialNum)',ones(trialNum,1); (1:trialNum)',ones(trialNum,1)*2];
trialOrder = randperm(size(tempIndex,1));
trialIndex = tempIndex(trialOrder,:);

% set keyboard
KbName('UnifyKeyNames');
skipKey = KbName('space');
escape = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
rKey = KbName('r');

pageUp = KbName('pageup'); % increase binocular deviation
pageDown = KbName('pagedown'); % decrease binocular deviation

% initial
calF(deviation);

if testMode 
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference', 'SkipSyncTests', 0);
end

AssertOpenGL;
InitializeMatlabOpenGL;

screenId = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');

% Define background color:
blackBackground = BlackIndex(screenId);

% Open a double-buffered full-screen window on the main displays screen.
[win , winRect] = PsychImaging('OpenWindow', screenId, blackBackground);
SCREEN.widthPix = winRect(3);
SCREEN.heightPix = winRect(4);

SCREEN.center = [SCREEN.widthPix,SCREEN.heightPix]/2;
SCREEN.refreshRate = Screen('NominalFrameRate', screenId);
fP = [SCREEN.widthPix/2,SCREEN.heightPix/2];
fixation = [fP(1)-2 fP(2)-2 fP(1)+2 fP(2)+2];

Screen('BeginOpenGL', win);
glViewport(0, 0, RectWidth(winRect), RectHeight(winRect));
glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
Screen('EndOpenGL', win);

px = 0; py = 0; pz = SCREEN.distance;
vf = [0;0;-5];

GenerateStarField()

Screen('BlendFunction', win, GL_ONE, GL_ZERO);
Screen('FillRect', win ,blackBackground,[0 0 SCREEN.widthPix SCREEN.heightPix]);
Screen('Flip', win,0,0);
triali = 1;
choice = zeros(trialNum*2,1);
choiceTime = zeros(trialNum*2,1);
conditionIndex = nan(length(trialOrder),size(conditions,2)+1);

while triali <= size(trialIndex,1)
    trialCondition = conditions(trialIndex((triali),1),:);
    motionType = trialIndex(triali,2);
    
    [p,f,theta] = calMovement(trialCondition,motionType);
    
    for framei = 1:length(theta)
        [ ~, ~, keyCode] = KbCheck;
        if keyCode(escape)
            break;
        end
        
        if keyCode(pageUp)
            deviation = deviation + deviationAdjust;
            disp(['binocular deviation: ' num2str(deviation)]);
        elseif keyCode(pageDown)
            if deviation > deviationAdjust
                deviation = deviation - deviationAdjust;
                disp(['binocular deviation: ' num2str(deviation)]);
            else 
                deviation = 0;
                disp(['binocular deviation: ' num2str(deviation)]);
            end
        end
        
        calF(deviation,theta(framei));
        
        %% draw for left eye
        Screen('BeginOpenGL', win);
        glColorMask(GL.TRUE, GL.FALSE, GL.FALSE, GL.FALSE);
        glMatrixMode(GL.PROJECTION);
        glLoadIdentity;
        glFrustum( SCREEN.sinisterLeft,SCREEN.sinisterRight, SCREEN.bottom, SCREEN.top, SCREEN.near, SCREEN.far);
        glMatrixMode(GL.MODELVIEW);
        glLoadIdentity;
        gluLookAt(p(1,framei)-deviation/2,p(2,framei),p(3,framei),p(1,framei)-deviation/2+f(1,framei),p(2,framei)+f(2,framei),p(3,framei)+f(3,framei),0.0,1.0,0.0);
        glClearColor(0,0,0,0);
        glColor3f(1,1,0);
        
        % draw the fixation point and 3d dots
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        
        %% draw for right eye
        glColorMask(GL.FALSE, GL.TRUE, GL.FALSE, GL.FALSE);
        glMatrixMode(GL.PROJECTION);
        glLoadIdentity;
        glFrustum( SCREEN.dexterLeft,SCREEN.dexterRight, SCREEN.bottom, SCREEN.top, SCREEN.near, SCREEN.far);
        glMatrixMode(GL.MODELVIEW);
        glLoadIdentity;
        gluLookAt(p(1,framei)+deviation/2,p(2,framei),p(3,framei),p(1,framei)+deviation/2+f(1,framei),p(2,framei)+f(2,framei),p(3,framei)+f(3,framei),0.0,1.0,0.0);
        glClearColor(0,0,0,0);
        glColor3f(1,1,0);
        
        % draw the fixation point and 3d dots
        DrawDots3D(win,[STARDATA.x ; STARDATA.y; STARDATA.z]);
        
        Screen('EndOpenGL', win);
        
        Screen('FillOval', win, [255 255 0 255], fixation);
        
        Screen('Flip', win, 0, 0);
    end
    [ ~, ~, keyCode] = KbCheck;
    if keyCode(escape)
        break;
    end
    
    [~, ~, ~] = DrawFormattedText(win, 'What''s your heading direction?','center',SCREEN.center(2)/2,[200 200 200]);
    Screen('TextBackgroundColor',win, [0 0 0 0]);
    Screen('DrawingFinished',win);
    Screen('Flip',win,0,0);
    
    choiceT = tic;
    while toc(choiceT) <= choicePeriod
        [ ~, ~, keyCode] = KbCheck;
       if keyCode(leftKey)
           choice(triali) = 1;
           choiceTime(triali) = toc(choiceT);
           break
       elseif keyCode(rightKey)
           choice(triali) = 2;
           choiceTime(triali) = toc(choiceT);
           break
       end
    end
    if feedback
        degi = trialCondition(1);
        correctChoice = (degi >= 0)+1;
        if  choice(triali) == correctChoice
            sound(sin(2*pi*25*(1:3000)/200)); % correct cue
            [~, ~, ~] = DrawFormattedText(win, 'You are right!','center',SCREEN.center(2)/2,[20 200 20]);
        elseif choice(triali)
            sound(sin(2*pi*25*(1:3000)/600)); % wrong cue
            [~, ~, ~] = DrawFormattedText(win, 'Please try again.','center',SCREEN.center(2)/2,[200 20 20]);
        else
            sound(sin(2*pi*25*(1:3000)/600)); % missing cue
            [~, ~, ~] = DrawFormattedText(win, 'Oops, you missed this trial.','center',SCREEN.center(2)/2,[200 20 20]);
        end
        Screen('TextBackgroundColor',win, [0 0 0 0]);
        Screen('DrawingFinished',win);
        Screen('Flip',win,0,0);
        WaitSecs(feedbackDuration);
    end
    if choice(triali)
        conditionIndex(triali,:) = [trialCondition motionType];
        triali = triali+1;
    else
        trialIndex = cat(1,trialIndex,trialIndex(triali,:));
        trialIndex(triali,:) = [];
    end
end
save(fullfile(saveDir,fileName),'MOVE','conditions','trialIndex','conditionIndex','choice','choiceTime');
Screen('CloseAll');
cd(curdir);

function calF(deviation,deltaDegree)
global SCREEN
% calculate frustum
SCREEN.near = 20;
SCREEN.far = 100;
SCREEN.top = (SCREEN.near / SCREEN.distance) * (SCREEN.heightCM / 2.0);
SCREEN.bottom = (SCREEN.near / SCREEN.distance) * (-SCREEN.heightCM / 2.0);

if nargin ==1
    % left eye
    SCREEN.sinisterRight = (SCREEN.near / SCREEN.distance) * (SCREEN.widthCM / 2.0 + deviation / 2.0);
    SCREEN.sinisterLeft = (SCREEN.near / SCREEN.distance) * (-SCREEN.widthCM / 2.0 + deviation / 2.0);
    % right eye
    SCREEN.dexterRight = (SCREEN.near / SCREEN.distance) * (SCREEN.widthCM / 2.0 - deviation / 2.0);
    SCREEN.dexterLeft = (SCREEN.near / SCREEN.distance) * (-SCREEN.widthCM / 2.0 - deviation / 2.0);
elseif nargin == 2
    delta = SCREEN.near * sind(deltaDegree);
    % left eye
    SCREEN.sinisterRight = (SCREEN.near / SCREEN.distance) * (SCREEN.widthCM / 2.0 + deviation / 2.0)+delta;
    SCREEN.sinisterLeft = (SCREEN.near / SCREEN.distance) * (-SCREEN.widthCM / 2.0 + deviation / 2.0)+delta;
    % right eye
    SCREEN.dexterRight = (SCREEN.near / SCREEN.distance) * (SCREEN.widthCM / 2.0 - deviation / 2.0)+delta;
    SCREEN.dexterLeft = (SCREEN.near / SCREEN.distance) * (-SCREEN.widthCM / 2.0 - deviation / 2.0)+delta;
end
end

function DrawDots3D(windowPtr, xyz)
% Draw a large number of dots in 3D very efficiently.
%
% Usage: moglDrawDots3D(windowPtr, xyz);
%
% This function is the 3D equivalent of the Screen('DrawDots') subfunction
% for fast drawing of 2D dots. It has mostly the same paramters as that
% function, but extended into the 3D domain. It accepts a subset of the
% parameters for that function, ie., it is less liberal in what it accepts,
% in order to allow for simpler code and maximum efficiency.
%
% As a bonus, it accepts one additional parameter 'glslshader', the
% optional handle to a GLSL shading program, e.g., as loaded from
% filesystem via LoadGLSLProgram().
%
% The routine will draw into the 3D OpenGL userspace rendering context of
% onscreen window or offscreen window (or texture) 'windowPtr'. It will
% automatically switch to that window if it isn't already active in 3D
% mode, and it will restore the drawing target to whatever was set before
% invocation in whatever mode (2D or 3D). This is a convenience feature for
% lazy users that mostly draw in 2D. If you intend to draw more stuff in 3D
% for a given frame, then you should switch your targetwindow 'windowPtr'
% into 3D mode manually via Screen('BeginOpenGL') yourself beforehand. This
% will avoid redundant and expensive context switches and increase the
% execution speed of your code!
%
% Parameters and their meaning:
%
% 'windowPtr' Handle of window or texture to draw into.
% 'xyz' A 3-by-n or 4-by-n matrix of n dots to draw. Each column defines
% one dot to draw, either as 3D position (x,y,z) or 4D position (x,y,z,w).
% Must be a double matrix!
%
%
%
% 'dot_type' optional: A setting of zero will draw rectangular dots, a
% setting of 1 will draw round dots, a setting of 2 will draw round dots of
% extra high quality if the hardware supports that. For anti-aliased dots
% you must select a setting of 1 or 2 and enable alpha blending as well.
%
% 'glslshader' optional: If omitted, shading state is not changed. If set
% to zero, then the standard fixed function OpenGL pipeline is used, like
% in Screen('DrawDots') (under most circumstances). If a positive
% glslshader handle to a GLSL shading program object is provided, that
% shader program will be used. You can use this, e.g., to bind a custom vertex
% shader to perform complex per-dot calculations very fast on the GPU.
%
% See

% History:
% 03/01/2009  mk  Written.

% Need global GL definitions:
global GL;

% Child protection:
if isempty(GL)
    error('Need OpenGL mode to be enabled, which is not the case! Call InitializeMatlabOpenGL at start of your script first!');
end

if nargin < 1 || isempty(windowPtr)
    error('"windowPtr" window handle missing! This is required!');
end

if nargin < 2 || isempty(xyz)
    % xyz dot position matrix is empty! Nothing to do for us:
    return;
end

% Want single xyz vector as a 3 or 4 row, 1 column vector:
if (size(xyz, 1) == 1)
    % if (size(xyz, 1) == 1) && (ismember(size(xyz, 2), [3,4]))
    xyz = transpose(xyz);
end

nvc = size(xyz, 1);
nrdots = size(xyz, 2);

% % if ~(nvc == 3 || nvc == 4) || nrdots < 1
% %     error('"xyz" argument must have 3 or 4 rows for x,y,z or x,y,z,w components and at least 1 column for at least one dot to draw!');
% % end

% Is the OpenGL userspace context for this 'windowPtr' active, as required?
[previouswin, IsOpenGLRendering] = Screen('GetOpenGLDrawMode');
PreIsOpenGLRendering = IsOpenGLRendering;

% Our target window windowPtr already active?
if previouswin ~= windowPtr
    % No. Wrong window. OpenGL rendering for this window active?
    if IsOpenGLRendering
        % Yes. We need to disable OpenGL mode for that other window and
        % switch to our window:
        Screen('EndOpenGL', windowPtr);
        
        % Our window is attached, but it is in 2D mode, not 3D mode yet:
        IsOpenGLRendering = 0;
    end
end

% Either a different window than ours is bound in 2D mode, then OpenGL
% rendering is not yet active and we need to switch to our window and to
% OpenGL rendering.
%
% Or we just switched from a different window in 3D mode to our window in
% 2D mode. Then we need to switch our window into 3D mode.
%
% In both cases, IsOpenGLRendering == false will indicate this.
%
% A third option is that our wanted window is already active and 3D OpenGL
% mode is already active. In that case IsOpenGLRendering == true and we
% don't need to do anything to switch modes:
if ~IsOpenGLRendering
    % Either some window, or our window bound in 2D mode. Need to switch to
    % our window in 3D mode:
    Screen('BeginOpenGL', windowPtr);
end

% Ok our target window and userspace OpenGL rendering context is bound, we
% can setup and execute the actual drawing:

% Reset dot size to 1.0:
glPointSize(10);

% glColor3f(1,1,0);
glEnable(GL.POINT_SMOOTH);
% % % totalnum=size(xyz,2);
% % % for i=1:totalnum
% % %     glBegin(GL.TRIANGLES)
% % %     glVertex3d(0,0);
% % %     glVertex3d(256,0);
% % %     glVertex3d(128,256);
% % %     glEnd;
% % % end

% Enable fast rendering of arrays:
glEnableClientState(GL.VERTEX_ARRAY);

% Pass a pointer to the start of the point-coordinate array:
glVertexPointer(nvc, GL.DOUBLE, 0, xyz);
% glDrawArrays(GL.POINTS, 0, nrdots);
glDrawArrays(GL.TRIANGLES, 0, nrdots);
% Disable fast rendering of arrays:
glDisableClientState(GL.VERTEX_ARRAY);
glVertexPointer(nvc, GL.DOUBLE, 0, 0);
glDisable(GL.POINT_SMOOTH);

% Our work is done. If a different window than our target window was
% active, we'll switch back to that window and its state:
if previouswin ~= windowPtr
    % Different window was active before our invocation. Need to disable
    % our 3D mode and switch back to that window (in 2D mode):
    Screen('EndOpenGL', previouswin);
    
    % Was that window in 3D mode, i.e., OpenGL rendering for that window was active?
    if PreIsOpenGLRendering
        % Yes. We need to switch that window back into 3D OpenGL mode:
        Screen('BeginOpenGL', previouswin);
    end
else
    % Our window was active beforehand. Was it in 2D mode? In that case we
    % need to switch our window back to 2D mode. Otherwise we'll just stay
    % in 3D mode:
    if ~PreIsOpenGLRendering
        % Was in 2D mode. We need to switch back to 2D:
        Screen('EndOpenGL', windowPtr);
    end
end

% Switchback complete. The graphics system is the same state as it was
% before our invocation.
return;
end

function GenerateStarField()
global STARFIELD;
global STARDATA;
totalDots = round(STARFIELD.dimensionX*STARFIELD.dimensionY*STARFIELD.dimensionZ*STARFIELD.density);
baseX=rand(1,totalDots)*STARFIELD.dimensionX-STARFIELD.dimensionX/2.0;
baseY=rand(1,totalDots)*STARFIELD.dimensionY-STARFIELD.dimensionY/2.0;
% baseZ=rand(1,totalDots)*STARFIELD.dimensionZ-STARFIELD.dimensionZ/4.0;
baseZ=rand(1,totalDots)*STARFIELD.dimensionZ-STARFIELD.dimensionZ;
STARDATA.x=zeros(1,3*totalDots);
STARDATA.y=zeros(1,3*totalDots);
STARDATA.z=zeros(1,3*totalDots);
starSize = degree2length(STARFIELD.starSize);

%Vertex1
Vertex1X=baseX - starSize/2.0;
Vertex1Y=baseY - starSize/2.0;
Vertex1Z=baseZ;

%Vertex2
Vertex2X=baseX;
Vertex2Y=baseY + starSize/2.0;
Vertex2Z=baseZ;

%Vertex3
Vertex3X=baseX + starSize/2.0;
Vertex3Y=baseY - starSize/2.0;
Vertex3Z=baseZ;
j=1;

for i=1:totalDots
    STARDATA.x(j)=Vertex1X(i);
    STARDATA.y(j)=Vertex1Y(i);
    STARDATA.z(j)=Vertex1Z(i);
    j=j+1;
    STARDATA.x(j)=Vertex2X(i);
    STARDATA.y(j)=Vertex2Y(i);
    STARDATA.z(j)=Vertex2Z(i);
    j=j+1;
    STARDATA.x(j)=Vertex3X(i);
    STARDATA.y(j)=Vertex3Y(i);
    STARDATA.z(j)=Vertex3Z(i);
    j=j+1;
end
end

function length = degree2length(degree)
% this function convert degree value to pixel value
% it is better been used to calculte the pixel length from the central point
% On X axis: dir = 1; On Y axis: dir = 2
global SCREEN

a=SCREEN.widthPix / SCREEN.widthCM;
b=SCREEN.heightPix / SCREEN.heightCM;

if nargin == 1
    if abs(a-b)/min(a,b) < 0.05
        length = tand(degree) * SCREEN.distance;
    else
        error('Error in screen parameter or screen config, or you should define it is horizontal(1) / vertical(2).')
    end
end
end

function conditions = calculateCondition()
% [degree velcity duration rotationDegree initialDegree]
global MOVE

conRsId = [sort(repmat(MOVE.rotationDegree',length(MOVE.initialDegree),1),1),...
    repmat(MOVE.initialDegree',length(MOVE.rotationDegree),1)];
conDRsId = [sort(repmat(MOVE.duration',size(conRsId,1),1),1),...
    repmat(conRsId,length(MOVE.duration),1)];
conVDRsId = [sort(repmat(MOVE.velocity',size(conDRsId,1),1),1),...
    repmat(conDRsId,length(MOVE.velocity),1)];
conditions = [sort(repmat(MOVE.degree',size(conVDRsId,1),1),1),...
    repmat(conVDRsId,length(MOVE.degree),1)];
end

function [p,f,theta] = calMovement(trialCondition,motionType)
global SCREEN
% [degree velcity duration rotationDegree initialDegree]
degi = trialCondition(1);
veli = trialCondition(2);
duri = trialCondition(3);
rdi = trialCondition(4);
idi = trialCondition(5);

distance = veli*duri;
frameNum = duri*SCREEN.refreshRate;

p = [sind(degi)*distance/frameNum*(1:frameNum); zeros(1,frameNum);-cosd(degi)*distance/frameNum*(1:frameNum)];
if motionType == 1
    % rotation
    f = [sind(-abs(rdi)/rdi*idi);0;-cosd(-abs(rdi)/rdi*idi)];
    metrix = roty(-rdi/frameNum);
    for i = 2:frameNum
        f(:,i) = metrix*f(:,i-1);
    end
    theta = zeros(frameNum,1);
elseif motionType == 2
    % shift frustum
    f = repmat([0;0;-1],1,frameNum);
    theta  = rdi/frameNum*(1:frameNum)+idi;
end
end
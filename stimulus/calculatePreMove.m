function [pX,pY,pZ,fX,fY,fZ,thetaY] = calculatePreMove(motionType,index)
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.frameRate * TRIALINFO.preMoveDuration);

frameOrder = 1:frameNum;
pX = zeros(1,frameNum);
pY = zeros(1,frameNum);
pZ = (-TRIALINFO.headingSpeed .* TRIALINFO.preMoveDuration ./ frameNum) * frameOrder;

% extract facing direction
% if motionType ~= 3
    fX = pX;
    fY = pY;
    fZ = -SCREEN.distance * ones(1,frameNum);
    thetaY = zeros(1,frameNum);
% else
%     rotDeg = TRIALINFO.trialConditions{index(1)}(index(2),4);
%     fX = -sind(rotDeg/2)*SCREEN.distance*ones(1,frameNum);
%     fY = pY;
%     fZ = -cosd(rotDeg/2)*SCREEN.distance*ones(1,frameNum);
%     thetaY = ones(1,frameNum) * (-rotDeg/2);
% end
% end
    
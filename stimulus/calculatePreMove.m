function [pX,pY,pZ,fX,fY,fZ] = calculatePreMove()
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.refreshRate * TRIALINFO.preMoveDuration);

frameOrder = 1:frameNum;
pX = zeros(1,frameNum);
pY = zeros(1,frameNum);
pZ = (-TRIALINFO.headingSpeed .* TRIALINFO.preMoveDuration ./ frameNum) * frameOrder;


% extract facing direction
fX = pX;
fY = pY;
fZ = pZ - SCREEN.distance;
end
    
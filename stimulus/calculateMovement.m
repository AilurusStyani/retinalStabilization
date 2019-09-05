function [pX,pY,pZ,fX,fY,fZ,thetaY] = calculateMovement(motionType,index,pglX,pglY,pglZ,pfX,pfY,pfZ)
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.refreshRate * TRIALINFO.moveDuration);
headingDegree = TRIALINFO.trialConditions{index(1)}(index(2),2);
frameOrder = 1:frameNum;

% calculate position
pX = (TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum .* sind(headingDegree)) * frameOrder + pglX;
pY = zeros(1,frameNum) + pglY;
pZ =-(TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum .* cosd(headingDegree)) * frameOrder + pglZ;

if motionType == 3
    rotDeg = TRIALINFO.trialConditions{index(1)}(index(2),4);
    alpha =  rotDeg / frameNum;
    thetaY = alpha*(1:frameNum) - rotDeg/2;
else
    thetaY = 0*(1:frameNum);
end
fX = pfX*ones(1,frameNum);
fY = pfY*ones(1,frameNum);
fZ = pfZ*ones(1,frameNum);
end
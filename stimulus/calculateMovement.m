function [pX,pY,pZ,fX,fY,fZ,thetaY] = calculateMovement(motionType,index,pglX,pglY,pglZ)
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.refreshRate * TRIALINFO.moveDuration);
headingDegree = TRIALINFO.trialConditions{index(1)}(index(2),2);
frameOrder = 1:frameNum;

% calculate position
pX = (TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum .* sind(headingDegree)) * frameOrder + pglX;
pY = zeros(1,frameNum) + pglY;
pZ =-(TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum .* cosd(headingDegree)) * frameOrder + pglZ;


% calculate facing direction, look at the point on screen plane
pfX = ones(1,frameNum)*0;
pfY = ones(1,frameNum)*0;
pfZ = ones(1,frameNum)*(-SCREEN.distance);

if motionType == 3
    alpha = TRIALINFO.trialConditions{index(1)}(index(2),4) / frameNum;
    metrix = roty(-alpha);
    fMetrix = [pfX;pfY;pfZ];
    thetaY = alpha*(1:frameNum);
    for i = 2:frameNum
        fMetrix(:,i) = metrix * fMetrix(:,i-1);
    end
    pfX = fMetrix(1,:);
    pfY = fMetrix(2,:);
    pfZ = fMetrix(3,:);
else
    thetaY = 0*(1:frameNum);
end

fX = pfX+pX;
fY = pfY+pY;
fZ = pfZ+pZ;
end
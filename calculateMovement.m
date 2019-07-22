function [pX,pY,pZ,fX,fY,fZ] = calculateMovement(motionType,index,pglX,pglY,pglZ)
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.refreshRate * TRIALINFO.moveDuration);
headingDegree = TRIALINFO.trialConditions{index(1)}(index(2),2);
frameOrder = 1:frameNum;

% calculate position
pX = (TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum .* sind(headingDegree)) * frameOrder' + pglX;
pY = zeros(frameNum,1) + pglY;
pZ =-(TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum .* cosd(headingDegree)) * frameOrder' + pglZ;


% calculate facing direction
fX = pX;
fY = pY;
fZ = pZ-1;

if motionType == 3
    alpha = TRIALINFO.trialConditions{index(1)}(index(2),4) / frameNum;
    metrix = roty(-alpha);
    fMetrix = [(fX-pX)';(fY-pY)';(fZ-pZ)'];
    for i = 2:frameNum
        fMetrix(:,i) = metrix * fMetrix(:,i-1);
    end
    fX = fMetrix(1,:)'+pX;
    fY = fMetrix(2,:)'+pY;
    fZ = fMetrix(3,:)'+pZ;
end
end
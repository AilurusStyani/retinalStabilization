function [pX,pY,pZ,fX,fY,fZ] = calculateMovement(motionType,index,pglX,pglY,pglZ)
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.refreshRate * TRIALINFO.moveDuration);

frameOrder = 1:frameNum+1;
glX = zeros(frameNum+1,1);
glY = zeros(frameNum+1,1);
glZ = (-TRIALINFO.headingSpeed .* TRIALINFO.moveDuration ./ frameNum) * frameOrder';

% ectract position
pX = glX(1:end-1)+pglX;
pY = glY(1:end-1)+pglY;
pZ = glZ(1:end-1)+pglZ;

% extract facing direction
fX = glX(2:end);
fY = glY(2:end);
fZ = glZ(2:end);

if motionType == 3
    alpha = TRIALINFO.trialConditions{index(1)}(index(2),4) / frameNum;
    metrix = roty(-alpha);
    fMetrix = [fX';fY';fZ'];
    for i = 2:frameNum
        fMetrix(:,i) = metrix * fMetrix(:,i-1);
    end
    fX = fMetrix(1,:)'+pX;
    fY = fMetrix(2,:)'+pY;
    fZ = fMetrix(3,:)'+pZ;
else
    fX = fX+pglX;
    fY = fY+pglY;
    fZ = fZ+pglZ;
end
end
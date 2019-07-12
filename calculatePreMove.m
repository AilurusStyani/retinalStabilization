function [pX,pY,pZ,fX,fY,fZ] = calculatePreMove()
global TRIALINFO
global SCREEN

frameNum = ceil(SCREEN.refreshRate * TRIALINFO.preMoveDuration);

frameOrder = 1:frameNum+1;
glX = zeros(frameNum+1,1);
glY = zeros(frameNum+1,1);
glZ = (-TRIALINFO.headingSpeed .* TRIALINFO.preMoveDuration ./ frameNum) * frameOrder;

% ectract position
pX = glX(1:end-1);
pY = glY(1:end-1);
pZ = glZ(1:end-1);

% extract facing direction
fX = glX(2:end);
fY = glY(2:end);
fZ = glZ(2:end);
end
    
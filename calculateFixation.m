function [fX,fY] = calculateFixation()
global TRIALINFO
global SCREEN

SCREEN.widthPix

frameNum = ceil(SCREEN.refreshRate*TRIALINFO.moveDuration);
frameOrder = 1:frameNum;
fixPositionDX = -TRIALINFO.fixationInitialDegree + (TRIALINFO.fixSpeed*TRIALINFO.moveDuration/frameNum.*frameOrder);
fixPositionPX = degree2pix(fixPositionDX) + SCREEN.widthPix/2;
fixPositionPY = ones(frameNum,1)*SCREEN.heightPix/2;

fX{1} = fixPositionPX;
fX{2} = ones(frameNum,1)*SCREEN.widthPix/2;
fX{3} = -fixPositionPX;
fY{1} = fixPositionPY;
fY{2} = fixPositionPY;
fY{3} = fixPositionPY;
end
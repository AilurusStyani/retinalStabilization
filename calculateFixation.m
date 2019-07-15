function [fX,fY] = calculateFixation()
% fX{1} and fY{1} position of every frame for condition 1
% fX{2} and fY{2} position of every frame for condition 2
% fX{3} and fY{3} position of every frame for condition 3
global TRIALINFO
global SCREEN

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
function [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ)
% global SCREEN
global TRIALINFO

pXl = -TRIALINFO.deviation/2 + glX;
pYl = glY;
pZl = glZ;

fXl = fX; fYl = fY; fZl = fZ;

pXr = +TRIALINFO.deviation/2 + glX;
pYr = glY;
pZr = glZ;

fXr = fX; fYr = fY; fZr = fZ;
end

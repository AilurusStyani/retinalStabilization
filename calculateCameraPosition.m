function [pXl,pYl,pZl,fXl,fYl,fZl,pXr,pYr,pZr,fXr,fYr,fZr] = calculateCameraPosition(glX,glY,glZ,fX,fY,fZ)
global SCREEN
global TRIALINFO

alpha = atand((TRIALINFO.deviation/2)/SCREEN.distance);
metrixRY = roty(alpha);
metrixLY = roty(-alpha);
vectorf2p = [glX;glY;glZ] - [fX;fY;fZ]; % the vector from face direction (the point on screen plane) to camera position

% calculate for left eye
vectorf2pL = metrixLY *  vectorf2p;

pXl = vectorf2pL(1,:) + fX;
pYl = vectorf2pL(2,:) + fY;
pZl = vectorf2pL(3,:) + fZ;

fXl = fX;
fYl = fY;
fZl = fZ;

% calculate for right eye
vectorf2pR = metrixRY * vectorf2p;

pXr = vectorf2pR(1,:) + fX;
pYr = vectorf2pR(2,:) + fY;
pZr = vectorf2pR(3,:) + fZ;

fXr = fX;
fYr = fY;
fZr = fZ;
end

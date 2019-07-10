function [pX,pY,pZ] =calculateLinearMovementForFix(distance,time,towards)

global TRIALINFO;

v = distance/time;
step = 1.0/TRIALINFO.refreshRate;
frameNum = TRIALINFO.refreshRate*time + 1;
pX = zeros(frameNum,1);
pY = zeros(frameNum,1);
pZ = zeros(frameNum,1);

for i=1:frameNum-1
    
   disi = i*step*v;
   
   if(towards>0)
       pX(i+1) = disi;
   else
       pX(i+1) = -disi;
   end
       
end
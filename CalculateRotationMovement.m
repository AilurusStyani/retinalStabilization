function [angleTrajectory,rotVector]= CalculateRotationMovement(elevation,azimuth,rdegree,time)

 global TRIALINFO;   
DEG2RAD = 0.017453292519943;
frameNum= TRIALINFO.refreshRate*time+1;

rotationFrame= rdegree/frameNum;

angleTrajectory=zeros(frameNum,1);

elevation = elevation*DEG2RAD;
azimuth = azimuth*DEG2RAD;
cosE = cos(elevation);
cosA = cos(azimuth);
sinE = sin(elevation);
sinA = sin(azimuth);
rotVector(1)=cosE*cosA;
rotVector(2)=sinE;
rotVector(3)=-cosE*sinA;
for i=2:frameNum
   angleTrajectory(i+1) = (i-1)*  rotationFrame;
end

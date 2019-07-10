function [glX,glY,glZ]= CalculateConstentVelocityMovement(elevation,azimuth,distance,time)

 global TRIALINFO;   
 deg2rad = 0.017453292519943;
v = (distance)/time;
step = 1.0/TRIALINFO.refreshRate;
frameNum = TRIALINFO.refreshRate*time + 1;
glX=zeros(frameNum,1);
glY=zeros(frameNum,1);
glZ=zeros(frameNum,1);
elevation = elevation*deg2rad;
azimuth = azimuth*deg2rad;
for i=1:frameNum-1
   disi = i*step*v;
   glZ(i+1) = disi * sin(elevation);
   glY(i+1) = disi * cos(elevation) * sin(azimuth);
   glX(i+1) = disi * cos(elevation) * cos(azimuth);
end



function [glX,glY,glZ]= CalculateConstentVelocityMovement(elevation,azimuth,distance,time,refreshRate)

v = (distance)/time;
step = 1.0/refreshRate;
frameNum = refreshRate*time + 1;
glX=zeros(frameNum,1);
glY=zeros(frameNum,1);
glZ=zeros(frameNum,1);
for i=1:frameNum-1
   disi = i*step*v;
   glZ(i+1) = disi * sind(elevation);
   glY(i+1) = disi * cosd(elevation) * sind(azimuth);
   glX(i+1) = disi * cosd(elevation) * cosd(azimuth);
end
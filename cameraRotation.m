function [pX,pY,pZ,rX,rY,rZ] = cameraRotation(trialCondition)
global CAMERA
global TRIALINFO;
duration = trialCondition(5);
heading_azimuth = -trialCondition(1);
distance = trialCondition(2);
rotationDegree = trialCondition(3);

frameNum = round(duration * TRIALINFO.refreshRate);
degreeFrame = rotationDegree / frameNum;
distanceFrame = distance / frameNum;

glX = zeros(frameNum+1,1);
glY = ones(frameNum+1,1)*CAMERA.elevation;
glZ = zeros(frameNum+1,1);

rX = zeros(frameNum,1);
rY = ones(frameNum,1)*CAMERA.elevation;
rZ = zeros(frameNum,1);

% rotationAngle = zeros(frameNum,1);

for i = 1:frameNum
        glX(i+1) = glX(i) + distanceFrame*cosd(heading_azimuth);
        glZ(i+1) = glZ(i) - distanceFrame*sind(heading_azimuth);
end

% calculate the camera position
pX = glX(1:end-1);
pY = glY(1:end-1);
pZ = glZ(1:end-1);

% % % calculate the face forward of the camera
fX = glX(2:end);
fY = glY(2:end);
fZ = glZ(2:end);
if degreeFrame>0
    rMatrix = roty(-degreeFrame);
elseif degreeFrame<0
    rMatrix = roty(degreeFrame);
end

% rFace = [pX(1) pY(1) pZ(1)]; %%% direction of the first frame
% % rFace = [fX fY fZ];%%%direction of the original frame
rFace = [fX fY fZ];
for i = 1: frameNum-1
    
% % %     if degreeFrame>0
% % %         rMatrix = roty(-degreeFrame*i);
% % %     elseif degreeFrame<0
% % %         rMatrix = roty(degreeFrame*i);
% % %     end
    
    rFacei = rMatrix * rFace(i,:)';
  
          rX(i+1) = rFacei(1) + fX(i) ;
          rZ(i+1) = rFacei(3) + fZ(i) ;
    
         
  
    rFace(i+1,:) = [rX(i+1) rY(i+1) rZ(i+1)];
end


end


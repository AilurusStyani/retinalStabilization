function generateStars(position,starSize)
global STARDATA
starNum = size(position,1);
STARDATA.x=zeros(1,3*starNum);
STARDATA.y=zeros(1,3*starNum);
STARDATA.z=zeros(1,3*starNum);

j=1;
for i = 1:starNum
    STARDATA.x(j)=position(i,1)-starSize/2;
    STARDATA.y(j)=position(i,2)-starSize/2;
    STARDATA.z(j)=position(i,3);
    j=j+1;
    STARDATA.x(j)=position(i,1);
    STARDATA.y(j)=position(i,2)+starSize/2;
    STARDATA.z(j)=position(i,3);
    j=j+1;
    STARDATA.x(j)=position(i,1)+starSize/2;
    STARDATA.y(j)=position(i,2)-starSize/2;
    STARDATA.z(j)=position(i,3);
    j=j+1;
end
end
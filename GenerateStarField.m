function GenerateStarField()
 global STARFIELD;
 global STARDATA;
totalDots = round(STARFIELD.dimensionX*STARFIELD.dimensionY*STARFIELD.dimensionZ*STARFIELD.density);
baseX=rand(1,totalDots)*dimensionX-dimensionX/2.0;
baseY=rand(1,totalDots)*dimensionY-dimensionY/2.0;
baseZ=rand(1,totalDots)*dimensionZ-dimensionZ/4.0;
STARDATA.x=zeros(1,3*totalDots);
STARDATA.y=zeros(1,3*totalDots);
STARDATA.z=zeros(1,3*totalDots);
size = degree2pix(STARFIELD.starSize);

%Vertex1
Vertex1X=baseX - size/2.0;
Vertex1Y=baseY - size/2.0;
Vertex1Z=baseZ;

%Vertex2
Vertex2X=baseX;
Vertex2Y=baseY + size/2.0;
Vertex2Z=baseZ;

%Vertex3
Vertex3X=baseX + size/2.0;
Vertex3Y=baseY - size/2.0;
Vertex3Z=baseZ;
j=1;

for i=1:totalDots
    STARDATA.x(j)=Vertex1X(i);
    STARDATA.y(j)=Vertex1Y(i);
    STARDATA.z(j)=Vertex1Z(i);
    j=j+1;
    STARDATA.x(j)=Vertex2X(i);
    STARDATA.y(j)=Vertex2Y(i);
    STARDATA.z(j)=Vertex2Z(i);
    j=j+1;
    STARDATA.x(j)=Vertex3X(i);
    STARDATA.y(j)=Vertex3Y(i);
    STARDATA.z(j)=Vertex3Z(i);
    j=j+1;
end


function pixel = degree2pix(degree,dir)
% this function convert degree value to pixel value
% it is better been used to calculte the pixel length from the fixation point
% On X axis: dir = 1; On Y axis: dir = 2
global SCREEN

a=SCREEN.widthPix / SCREEN.widthCM;
b=SCREEN.heightPix / SCREEN.heightCM;

if nargin == 1
    if abs(a-b)/min(a,b) < 0.05
        length = tand(degree) * SCREEN.distance;
        pixel = length / SCREEN.widthCM * SCREEN.widthPix;
    else
        error('Error in screen parameter or screen config.')
    end
elseif nargin == 2
    length = tand(degree) * SCREEN.distance;
    if dir == 1
        pixel = length / SCREEN.widthCM * SCREEN.widthPix;
    elseif dir == 2
        pixel = length / SCREEN.heightCM * SCREEN.heightPix;
    else
        error('Invalid value for dir.')
    end
end
function degree = pix2degree(pix,dir)
% this function convert degree value to pixel value
% it is better been used to calculte the pixel length from the central point
% On X axis: dir = 1; On Y axis: dir = 2
global SCREEN

a=SCREEN.widthPix / SCREEN.widthCM;
b=SCREEN.heightPix / SCREEN.heightCM;

if nargin == 1
    if abs(a-b)/min(a,b) < 0.05
        if length(pix)==1
            l = pix * SCREEN.widthCM / SCREEN.widthPix;
            degree = atand(l/SCREEN.distance);
        elseif length(pix)==2
            lx = pix(1) * SCREEN.widthCM / SCREEN.widthPix;
            degreex = atand(lx/SCREEN.distance);
            ly = pix(2) * SCREEN.heightCM / SCREEN.heightPix;
            degreey = atand(ly/SCREEN.distance);
            degree = [degreex,degreey];
        end
    else
        error('Error in screen parameter or screen config, or you should define it is horizontal(1) / vertical(2).')
    end
elseif nargin == 2
    if dir == 1 && length(pix)==1
        l = pix * SCREEN.widthCM / SCREEN.widthPix;
        degree = atand(l/SCREEN.distance);
    elseif dir == 2 && length(pix)==1
        l = pix * SCREEN.heightCM / SCREEN.heightPix;
        degree = atand(l/SCREEN.distance);
    else
        error('Invalid value for input.')
    end
end
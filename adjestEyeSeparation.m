function adjestEyeSeparation(keyPlus,keyReduce,adjustValue)
global CAMERA
global FRUSTUM
global SCREEN

if keyPlus
    CAMERA.eyeSeparation = CAMERA.eyeSeparation + adjustValue;
elseif keyReduce
    CAMERA.eyeSeparation = CAMERA.eyeSeparation - adjustValue;
end

% right eye
FRUSTUM.oculusDexterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-(SCREEN.widthCM+CAMERA.eyeSeparation) / 2.0 );
FRUSTUM.oculusDexterRight = (FRUSTUM.clipNear / SCREEN.distance) * ((SCREEN.widthCM-CAMERA.eyeSeparation) / 2.0 );

% left eye
FRUSTUM.oculusSinisterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-(SCREEN.widthCM-CAMERA.eyeSeparation) / 2.0 );
FRUSTUM.oculusSinisterRight = (FRUSTUM.clipNear / SCREEN.distance) * ((SCREEN.widthCM+CAMERA.eyeSeparation) / 2.0 );
end
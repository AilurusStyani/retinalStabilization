function calculateFrustum(coordinateMuilty)
global TRIALINFO;
global FRUSTUM;
global SCREEN;

FRUSTUM.clipNear = SCREEN.distance; % cm
FRUSTUM.clipFar = 150*coordinateMuilty; % cm
FRUSTUM.top = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.heightCM / 2.0);
FRUSTUM.bottom = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.heightCM / 2.0);

% left eye
FRUSTUM.sinisterRight = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 + TRIALINFO.deviation / 2.0);
FRUSTUM.sinisterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 + TRIALINFO.deviation / 2.0);

% right eye
FRUSTUM.dexterRight = (FRUSTUM.clipNear / SCREEN.distance) * (SCREEN.widthCM / 2.0 - TRIALINFO.deviation / 2.0);
FRUSTUM.dexterLeft = (FRUSTUM.clipNear / SCREEN.distance) * (-SCREEN.widthCM / 2.0 - TRIALINFO.deviation / 2.0);
end
function retryFlag = pursuitCheck(pursuitPoint,eyeTrackWinP)

evt = Eyelink( 'NewestFloatSample');
eyeUsed = Eyelink('EyeAvailable'); % get eye that's tracked
if eyeUsed ~= -1 % do we know which eye to use yet?
    px =evt.gx(eyeUsed+1); % +1 as we're accessing MATLAB array
    py = evt.gy(eyeUsed+1);
end

if abs((pursuitPoint(1)-px)+(pursuitPoint(2)-py)*1i) > eyeTrackWinP
    retryFlag=1;
else
    retryFlag=0;
end
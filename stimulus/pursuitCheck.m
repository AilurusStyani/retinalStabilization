function retryFlag = pursuitCheck(pursuitPoint,eyeTrackWinP,win)
global TRIALINFO

if Eyelink( 'NewFloatSampleAvailable')>0
    % get the sample in the form of an event structure
    evt = Eyelink( 'NewestFloatSample');
    eyeUsed = Eyelink('EyeAvailable'); % get eye that's tracked
    if eyeUsed ~= -1 % do we know which eye to use yet?
        px =evt.gx(eyeUsed+1); % +1 as we're accessing MATLAB array
        py = evt.gy(eyeUsed+1);
    end
%     drawFixation([px,py],TRIALINFO.fixationSizeP,win);
    if abs((pursuitPoint(1)-px)+(pursuitPoint(2)-py)*1i) > eyeTrackWinP
        retryFlag=1;
    else
        retryFlag=0;
    end
else
    error('no new eye data')
end
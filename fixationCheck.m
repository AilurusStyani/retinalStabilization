function [escFlag,retryFlag] = fixationCheck(fixationPoint,eyeTrackerWinP,fixationPeriod,eye_used,escape,skipKey,cKey,el)
% fixationPoint: screen position of fixation point
% eyeTrackerWinP: fixation window of pixel
% fixationPeriod: the time window of fixation
% eye_used: in case not using left eye, this function is coded for left eye currently
% escape: break the block
% skipKey: skip the fixation
% cKey: force calibration
% el: eyelink stract
escFlag = 0;
retryFlag = 0;
while 1
    fixationStart = tic;
    while 1
        [~,~,keyCode] = KbCheck;
        if keyCode(escape)
            escFlag = 1;
            return
        elseif keyCode(skipKey)
            return
        end
        
        if keyCode(cKey) % press C to calibrate
            EyelinkDoTrackerSetup(el);
            Eyelink('StartRecording');
            eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
            if eye_used == el.BINOCULAR % if both eyes are tracked
                eye_used = el.LEFTEYE; % use left eye
            end
            Eyelink('message', 'SYNCTIME');	 	 % zero-plot time for EDFVIEW
            error=Eyelink('checkrecording'); 		% Check recording status */
            if(error~=0)
                fprintf('Eyelink checked wrong status.\n');
                cleanup;  % cleanup function
                Eyelink('ShutDown');
                Screen('CloseAll');
            end
            WaitSecs(0.5); % wait a bit
            retryFlag = 1;
            return
        end
        
%         if Eyelink( 'NewFloatSampleAvailable')>0
            % get the sample in the form of an event structure
            evt = Eyelink( 'NewestFloatSample');
            if eye_used ~= -1 % do we know which eye to use yet?
                px =evt.gx(eye_used+1); % +1 as we're accessing MATLAB array
                py = evt.gy(eye_used+1);
                % frameStartTime(i) = evt.time;
            end
%         end
        
        if abs((fixationPoint(1)-px)+(fixationPoint(2)-py)*1i) > eyeTrackerWinP
            break
        elseif toc(fixationStart) >= fixationPeriod
            return
        end
    end
end
end
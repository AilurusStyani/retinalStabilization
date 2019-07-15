function drawFixation(fixationPosition)
global TRIALINFO

sizeP = degree2pix(TRIALINFO.fixzationSizeD/2);
fixation = [fixationPosition(1)-sizeP, fixationPosition(2)-sizeP, fixationPosition(1)+sizeP, fixationPosition(2)+sizeP];
Screen('FillOval', win, [255 0 0 255], fixation);
end
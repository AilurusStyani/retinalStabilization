close all;
clear all;

filePath = 'D:\ION\2019\2019Internship\Dhwani\TestData';
ascFile = dir(fullfile(filePath,'*.asc'));
matFile = strrep(ascFile.name,'.asc','.mat');
data = load(fullfile(filePath,matFile));

trialLength = length(data.Conditions);

% motion type 1:
headDeg{1} = nan(trialLength,2);
headSpe{1} = nan(trialLength,2);

% motion type 2:
headDeg{2} = nan(trialLength,2);
headSpe{2} = nan(trialLength,2);
fixDir{2} = nan(trialLength,2);
fixSpe{2} = nan(trialLength,2);

% motion type 3:
headDeg{3} = nan(trialLength,2);
headSpe{3} = nan(trialLength,2);
rotationDeg{3} = nan(trialLength,2);
rotationSpe{3} = nan(trialLength,2);

% motion type 4:
headDeg{4} = nan(trialLength,2);
headSpe{4} = nan(trialLength,2);
fixDir{4} = nan(trialLength,2);
fixSpe{4} = nan(trialLength,2);

for i = 1:trialLength
    
    if isempty(data.Conditions{i})
        continue
    else
        motionType = data.Conditions{i}(end);
    end
    
    if motionType == 1
        headDeg{1}(i,:) = [data.Conditions{i}(2),i];
        headSpe{1}(i,:) = [data.Conditions{i}(3),i];
    elseif motionType == 2
        headDeg{2}(i,:) = [data.Conditions{i}(2),i];
        headSpe{2}(i,:) = [data.Conditions{i}(3),i];
        fixDir{2}(i,:) = [data.Conditions{i}(4),i];
        fixSpe{2}(i,:) = [data.Conditions{i}(5),i];
    elseif motionType == 3
        headDeg{3}(i,:) = [data.Conditions{i}(2),i];
        headSpe{3}(i,:) = [data.Conditions{i}(3),i];
        rotationDeg{3}(i,:) = [data.Conditions{i}(4),i];
        rotationSpe{3}(i,:) = [data.Conditions{i}(5),i];
    elseif motionType == 4
        headDeg{4}(i,:) = [data.Conditions{i}(2),i];
        headSpe{4}(i,:) = [data.Conditions{i}(3),i];
        fixDir{4}(i,:) = [data.Conditions{i}(4),i];
        fixSpe{4}(i,:) = [data.Conditions{i}(5),i];
    end
end

uniqueDeg = cell(1,4);
for j = 1:4
    del = isnan(headDeg{j}(:,1));
    headDeg{j}(del,:) = [];
    del = isnan(headSpe{j}(:,1));
    headSpe{j}(del,:) = [];
    if j == 2 || j == 4
        del = isnan(fixDir{j}(:,1));
        fixDir{j}(del,:) = [];
        del = isnan(fixSpe{j}(:,1));
        fixSpe{j}(del,:) = [];
    elseif j ==3
        del = isnan(rotationDeg{j}(:,1));
        rotationDeg{j}(del,:) = [];
        del = isnan(rotationSpe{j}(:,1));
        rotationSpe{j}(del,:) = [];
    end
    uniqueDeg{j} = unique(headDeg{j}(:,1));
end

% [~,hSubplot1] = tight_subplot(1,4,[0.02 0.02]);
for j=1:4
    choice = cell(length(uniqueDeg{j}),1);
    pRight = cell(length(uniqueDeg{j}),1);
    choiceTime = cell(length(uniqueDeg{j}),1);
    for k = 1:length(uniqueDeg{j})
        index = ismember(headDeg{j}(:,1),uniqueDeg{j}(k));
        trialIndex = headDeg{j}(index,2);
        choice{k} = data.choice(trialIndex);
        pRight{k} = sum(choice{k} == 2)./ length(choice{k});
        choiceTime{k} = length(choice{k});
    end
    fitData = [uniqueDeg{j},cell2mat(pRight),cell2mat(choiceTime)];
    
    [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
    xi = min(uniqueDeg{j}):0.1:max(uniqueDeg{j});
    y_fit = cum_gaussfit([bias,threshold],xi);
    
    if ishandle(j); end; figure(j); set(gcf,'color','white');
    plot([90,90],[0,1],'-.k');
    hold on
    plot(uniqueDeg{j},cell2mat(pRight),'*');
    plot(xi,y_fit,'-');
    set(gca, 'xlim',[-50,50])
    xlabel('Heading degree');
    ylabel('Proportion of "right" choice');
%     hleg1=legend('choice','mean & standard error','linear result');
%     set(hleg1,'Location','EastOutside')
    text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
    text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
end

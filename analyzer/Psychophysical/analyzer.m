close all;
clear all;

filePath = 'D:\BYC\2019Intern\2019Internship\Dhwani\TestData';
edfFile = dir(fullfile(filePath,'*.edf'));
flipNameStr = 'flip';

popChoice = [];
for i = 1:length(edfFile)
    matFile = strrep(edfFile(i).name,'.edf','.mat');
    data = load(fullfile(filePath,matFile));
    popChoice = cat(1,popChoice,data.choice);
    trialLength = length(data.Conditions);
    flipDegree = contains(matFile,flipNameStr);
    
    % motion type 1:
    headDeg{1} = nan(trialLength,2);
    headSpe{1} = nan(trialLength,2);
    popHeadDeg{1} = [];
    popHeadSpe{1} = [];
    
    % motion type 2:
    headDeg{2} = nan(trialLength,2);
    headSpe{2} = nan(trialLength,2);
    fixDir{2} = nan(trialLength,2);
    fixSpe{2} = nan(trialLength,2);
    popHeadDeg{2} = [];
    popHeadSpe{2} = [];
    popFixDir{2} = [];
    popFixSpe{2} = [];
    
    % motion type 3:
    headDeg{3} = nan(trialLength,2);
    headSpe{3} = nan(trialLength,2);
    rotationDeg{3} = nan(trialLength,2);
    rotationSpe{3} = nan(trialLength,2);
    popHeadDeg{3} = [];
    popHeadSpe{3} = [];
    popRotationDeg{3} = [];
    popRotationSpe{3} = [];
    
    % motion type 4:
    headDeg{4} = nan(trialLength,2);
    headSpe{4} = nan(trialLength,2);
    fixDir{4} = nan(trialLength,2);
    fixSpe{4} = nan(trialLength,2);
    popHeadDeg{4} = [];
    popHeadSpe{4} = [];
    popFixDir{4} = [];
    popFixSpe{4} = [];
    
    for j = 1:trialLength
        
        if isempty(data.Conditions{j})
            continue
        else
            motionType = data.Conditions{j}(end);
        end
        
        if motionType == 1
            headDeg{1}(j,:) = [data.Conditions{j}(2),j];
            headSpe{1}(j,:) = [data.Conditions{j}(3),j];
        elseif motionType == 2
            headDeg{2}(j,:) = [data.Conditions{j}(2),j];
            headSpe{2}(j,:) = [data.Conditions{j}(3),j];
            fixDir{2}(j,:) = [data.Conditions{j}(4),j];
            fixSpe{2}(j,:) = [data.Conditions{j}(5),j];
        elseif motionType == 3
            headDeg{3}(j,:) = [data.Conditions{j}(2),j];
            headSpe{3}(j,:) = [data.Conditions{j}(3),j];
            rotationDeg{3}(j,:) = [data.Conditions{j}(4),j];
            rotationSpe{3}(j,:) = [data.Conditions{j}(5),j];
        elseif motionType == 4
            headDeg{4}(j,:) = [data.Conditions{j}(2),j];
            headSpe{4}(j,:) = [data.Conditions{j}(3),j];
            fixDir{4}(j,:) = [data.Conditions{j}(4),j];
            fixSpe{4}(j,:) = [data.Conditions{j}(5),j];
        end
    end
    
    uniqueDeg = cell(1,4);
    for j = 1:4
        del = isnan(headDeg{j}(:,1));
        headDeg{j}(del,:) = [];
        
        if flipDegree
            headDeg{j} = -headDeg{j};
        end
        
        del = isnan(headSpe{j}(:,1));
        headSpe{j}(del,:) = [];
        if j == 2 || j == 4
            del = isnan(fixDir{j}(:,1));
            fixDir{j}(del,:) = [];
            del = isnan(fixSpe{j}(:,1));
            fixSpe{j}(del,:) = [];
            
            popFixDir{j} = cat(1,popFixDir{j},fixDir{j});
            popFixSpe{j} = cat(1,popFixSpe{j},fixSpe{j});
        elseif j ==3
            del = isnan(rotationDeg{j}(:,1));
            rotationDeg{j}(del,:) = [];
            del = isnan(rotationSpe{j}(:,1));
            rotationSpe{j}(del,:) = [];
            
            popRotationDeg{j} = cat(1,popRotationDeg{j},rotationDeg{j});
            popRotationSpe{j} = cat(1,popRotationSpe{j},rotationSpe{j});
        end
        uniqueDeg{j} = unique(headDeg{j}(:,1));
        
        popHeadDeg{j} = cat(1,popHeadDeg{j},headDeg{j});
        popHeadSpe{j} = cat(1,popHeadSpe{j},headSpe{j});
    end
    
    % [~,hSubplot1] = tight_subplot(1,4,[0.02 0.02]);
    for j=1:4
        choice = cell(length(uniqueDeg{j}),1);
        pRight = cell(length(uniqueDeg{j}),1);
        choiceTime = cell(length(uniqueDeg{j}),1);
        for k = 1:length(uniqueDeg{j})
            index = ismember(headDeg{j}(:,1),uniqueDeg{j}(k));
            trialIndex = headDeg{j}(index,2);
            choice{k} = data.choice(trialIndex); %%
            pRight{k} = sum(choice{k} == 2)./ length(choice{k});
            choiceTime{k} = length(choice{k});
        end
        fitData = [uniqueDeg{j},cell2mat(pRight),cell2mat(choiceTime)];
        
        [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
        xi = min(uniqueDeg{j}):0.1:max(uniqueDeg{j});
        y_fit = cum_gaussfit([bias,threshold],xi);
        
        if ishandle((i-1)*length(edfFile)+j); end; figure((i-1)*length(edfFile)+j); set(gcf,'color','white');
        plot([90,90],[0,1],'-.k');
        hold on
        plot(uniqueDeg{j},cell2mat(pRight),'*');
        plot(xi,y_fit,'-');
        set(gca, 'xlim',[-15,15],'ylim',[0 1])
        xlabel('Heading degree');
        ylabel('Proportion of "right" choice');
        title(['Result for condition ' num2str(j)]);
        %     hleg1=legend('choice','mean & standard error','linear result');
        %     set(hleg1,'Location','EastOutside')
        text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
        text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
    end
end

% draw popular result
popUniqueDeg = cell(1,4);
for i = 1:4
    popUniqueDeg{i} = unique(popHeadDeg{i}(:,1));
    choice = cell(length(popUniqueDeg{i}),1);
    pRight = cell(length(popUniqueDeg{i}),1);
    choiceTime = cell(length(popUniqueDeg{i}),1);
    for k = 1:length(popUniqueDeg{i})
        index = ismember(popHeadDeg{i}(:,1),popUniqueDeg{i}(k));
        trialIndex = popHeadDeg{i}(index,2);
        choice{k} = data.choice(trialIndex);
        pRight{k} = sum(choice{k} == 2)./ length(choice{k});
        choiceTime{k} = length(choice{k});
    end
    fitData = [popUniqueDeg{i},cell2mat(pRight),cell2mat(choiceTime)];
    
    [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
    xi = min(popUniqueDeg{i}):0.1:max(popUniqueDeg{i});
    y_fit = cum_gaussfit([bias,threshold],xi);
    
    if ishandle(4*length(edfFile)+i); end; figure(4*length(edfFile)+i); set(gcf,'color','white');
    plot([90,90],[0,1],'-.k');
    hold on
    plot(popUniqueDeg{i},cell2mat(pRight),'*');
    plot(xi,y_fit,'-');
    set(gca, 'xlim',[-15,15],'ylim',[0 1]);
    title(['Popular result for condition ' num2str(i)]);
    xlabel('Heading degree');
    ylabel('Proportion of "right" choice');
    %     hleg1=legend('choice','mean & standard error','linear result');
    %     set(hleg1,'Location','EastOutside')
    text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
    text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
end


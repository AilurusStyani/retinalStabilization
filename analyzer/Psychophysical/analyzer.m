close all;
clear all;

filePath = 'D:\BYC\2019Intern\2019Internship\Dhwani\retinalStabilization\stimulus\data';
edfFile = dir(fullfile(filePath,'*.edf'));
flipNameStr = 'flip';
pBias = nan(4,length(edfFile));
pThreshold = nan(4,length(edfFile));
popChoice = [];
pI = 0;
% hPopulation = tight_subplot(2,2,[0.5 0.5]);
% set(hPopulation,'xtick',[]);
% set(hPopulation,'ytick',[]);

popHeadDeg{1} = [];
popHeadSpe{1} = [];

popHeadDeg{2} = [];
popHeadSpe{2} = [];
popPursuitDir{2} = [];
popFixSpe{2} = [];

popHeadDeg{3} = [];
popHeadSpe{3} = [];
popPursuitDir{3} = [];
popRotationSpe{3} = [];

popHeadDeg{4} = [];
popHeadSpe{4} = [];
popPursuitDir{4} = [];
popFixSpe{4} = [];

pChoice = cell(4,1);

figureNum = 1;
for i = 1:length(edfFile)
    matFile = strrep(edfFile(i).name,'.edf','.mat');
    data = load(fullfile(filePath,matFile));
    popChoice = cat(1,popChoice,data.choice);
    trialLength = length(data.Conditions);
    flipDegree = contains(matFile,flipNameStr);
    num = regexp(matFile, '(\d+)','tokens');
    subNum = cell2mat(num{1});
    
    % motion type 1:
    headDeg{1} = nan(trialLength,2);
    headSpe{1} = nan(trialLength,2);
    pursuitDir{1} = [];
    
    
    % motion type 2:
    headDeg{2} = nan(trialLength,2);
    headSpe{2} = nan(trialLength,2);
    pursuitDir{2} = nan(trialLength,2);
    fixSpe{2} = nan(trialLength,2);
    
    % motion type 3:
    headDeg{3} = nan(trialLength,2);
    headSpe{3} = nan(trialLength,2);
    pursuitDir{3} = nan(trialLength,2);
    rotationSpe{3} = nan(trialLength,2);
    
    % motion type 4:
    headDeg{4} = nan(trialLength,2);
    headSpe{4} = nan(trialLength,2);
    pursuitDir{4} = nan(trialLength,2);
    fixSpe{4} = nan(trialLength,2);
    
    
    
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
            pursuitDir{2}(j,:) = [data.Conditions{j}(4),j];
            fixSpe{2}(j,:) = [data.Conditions{j}(5),j];
        elseif motionType == 3
            headDeg{3}(j,:) = [data.Conditions{j}(2),j];
            headSpe{3}(j,:) = [data.Conditions{j}(3),j];
            pursuitDir{3}(j,:) = [data.Conditions{j}(4),j];
            rotationSpe{3}(j,:) = [data.Conditions{j}(5),j];
        elseif motionType == 4
            headDeg{4}(j,:) = [data.Conditions{j}(2),j];
            headSpe{4}(j,:) = [data.Conditions{j}(3),j];
            pursuitDir{4}(j,:) = [data.Conditions{j}(4),j];
            fixSpe{4}(j,:) = [data.Conditions{j}(5),j];
        end
    end
    
    uniqueDeg = cell(1,4);
    for j = 1:4
        del = isnan(headDeg{j}(:,1));
        headDeg{j}(del,:) = [];
        
        if flipDegree
            headDeg{j}(:,1) = -headDeg{j}(:,1);
        end
        
        del = isnan(headSpe{j}(:,1));
        headSpe{j}(del,:) = [];
        if j == 2 || j == 4
            del = isnan(pursuitDir{j}(:,1));
            pursuitDir{j}(del,:) = [];
            del = isnan(fixSpe{j}(:,1));
            fixSpe{j}(del,:) = [];
            
            popPursuitDir{j} = cat(1,popPursuitDir{j},pursuitDir{j});
            popFixSpe{j} = cat(1,popFixSpe{j},fixSpe{j});
        elseif j ==3
            del = isnan(pursuitDir{j}(:,1));
            pursuitDir{j}(del,:) = [];
            del = isnan(rotationSpe{j}(:,1));
            rotationSpe{j}(del,:) = [];
            
            popPursuitDir{j} = cat(1,popPursuitDir{j},pursuitDir{j});
            popRotationSpe{j} = cat(1,popRotationSpe{j},rotationSpe{j});
        end
        uniqueDeg{j} = unique(headDeg{j}(:,1));
        
        popHeadDeg{j} = cat(1,popHeadDeg{j},headDeg{j});
        popHeadSpe{j} = cat(1,popHeadSpe{j},headSpe{j});
        pChoice{j} = cat(1,pChoice{j},data.choice(headDeg{j}(:,2)));
    end
    
    % [~,hSubplot1] = tight_subplot(1,4,[0.02 0.02]);
    for j=1:4
        choice = cell(length(uniqueDeg{j}),1);
        pRight = cell(length(uniqueDeg{j}),1);
        choiceTime = cell(length(uniqueDeg{j}),1);
        
        if j > 1
            uniqueDir = unique(pursuitDir{j}(:,1));
            dirLength = length(uniqueDir);
            for dir = 1:dirLength
                trialIndexI = find(pursuitDir{j}(:,1) == uniqueDir(dir));
                dirIndex{dir} = pursuitDir{j}(trialIndexI,:);
            end
        end
        
        if j == 1
            for k = 1:length(uniqueDeg{j})
                index = ismember(headDeg{j}(:,1),uniqueDeg{j}(k));
                trialIndex = headDeg{j}(index,2);
                choice{k} = data.choice(trialIndex); %%
                pRight{k} = sum(choice{k} == 2)./ length(choice{k});
                choiceTime{k} = length(choice{k});
            end
            fitData = [uniqueDeg{j},cell2mat(pRight),cell2mat(choiceTime)];
            
            [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
            pBias(j,i) = bias;
            pThreshold(j,i) = threshold;
            xi = min(uniqueDeg{j}):0.1:max(uniqueDeg{j});
            y_fit = cum_gaussfit([bias,threshold],xi);
            
            if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
            plot([90,90],[0,1],'-.k');
            hold on
            plot(uniqueDeg{j},cell2mat(pRight),'*');
            plot(xi,y_fit,'-');
            set(gca, 'xlim',[-15,15],'ylim',[0 1])
            xlabel('Heading degree');
            ylabel('Proportion of "right" choice');
            title(['Condition ' num2str(j) ' for participant ' subNum]);
            %     hleg1=legend('choice','mean & standard error','linear result');
            %     set(hleg1,'Location','EastOutside')
            text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
            text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
            figureNum = figureNum+1;
        elseif j>1
            for dir = 1:dirLength
                for k = 1:length(uniqueDeg{j})
                    index = find(ismember(headDeg{j}(:,1),uniqueDeg{j}(k)));
                    trialIndex = headDeg{j}(index,2);
                    trialIndex = trialIndex(ismember(trialIndex,dirIndex{dir}(:,2)));
                    choice{k} = data.choice(trialIndex); %%
                    pRight{k} = sum(choice{k} == 2)./ length(choice{k});
                    choiceTime{k} = length(choice{k});
                    
                end
                fitData = [uniqueDeg{j},cell2mat(pRight),cell2mat(choiceTime)];
                if min(cell2mat(pRight)) == 1
                    fitData(:,2) = fitData(:,2) - 0.01;
                elseif max(cell2mat(pRight)) == 0
                    fitData(:,2) = fitData(:,2) + 0.01;
                end
                [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
                pBias(j,i) = bias;
                pThreshold(j,i) = threshold;
                xi = min(uniqueDeg{j}):0.1:max(uniqueDeg{j});
                y_fit = cum_gaussfit([bias,threshold],xi);
                
                if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
                plot([90,90],[0,1],'-.k');
                hold on
                plot(uniqueDeg{j},cell2mat(pRight),'*');
                plot(xi,y_fit,'-');
                set(gca, 'xlim',[-15,15],'ylim',[0 1])
                xlabel('Heading degree');
                ylabel('Proportion of "right" choice');
                title(['Condition ' num2str(j) ' direction ' num2str(dir) ' for participant ' subNum]);
                %     hleg1=legend('choice','mean & standard error','linear result');
                %     set(hleg1,'Location','EastOutside')
                text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
                text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
                figureNum = figureNum+1;
            end
        end
    end
end

% draw popular result
popUniqueDeg = cell(1,4);
colorIndex = [0 0 0; 1 0.0 0.0; 0 1 0; 0 0 1];
if ishandle(1000); end; figure(1000); clf; set(gcf,'color','white'); hold on;plot([0,0],[0,1],'-.k');
for i = 1:4
    popUniqueDeg{i} = unique(popHeadDeg{i}(:,1));
    choice = cell(length(popUniqueDeg{i}),1);
    pRight = cell(length(popUniqueDeg{i}),1);
    choiceTime = cell(length(popUniqueDeg{i}),1);
    pI = pI+length(popHeadDeg{i}(:,1));
    for k = 1:length(popUniqueDeg{i})
        index = find(popHeadDeg{i}(:,1) == popUniqueDeg{i}(k));
        choice{k} = pChoice{i}(index);
        pRight{k} = sum(choice{k} == 2)./ length(choice{k});
        choiceTime{k} = length(choice{k});
    end
    fitData = [popUniqueDeg{i},cell2mat(pRight),cell2mat(choiceTime)];
    
    [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
    xi = min(popUniqueDeg{i}):0.1:max(popUniqueDeg{i});
    y_fit = cum_gaussfit([bias,threshold],xi);
    
    
    hold on
    plot(popUniqueDeg{i},cell2mat(pRight),'*','color',colorIndex(i,:));
    h(i) = plot(xi,y_fit,'-','color',colorIndex(i,:));
    set(gca, 'xlim',[-15,15],'ylim',[0 1]);
    title('Population result');
    xlabel('Heading degree');
    ylabel('Proportion of "right" choice');
    %     hleg1=legend('choice','mean & standard error','linear result');
    %     set(hleg1,'Location','EastOutside')
    %     text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
    %     text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
end
legend(h,'fixation','normal pursuit','simulated pursuit','stabilized pursuit');
[h_Bias,p_Bias] = ttest2(pBias(1,:),pBias(3,:))
[h_Threshold,p_Threshold] = ttest2(pThreshold(1,:),pThreshold(3,:))

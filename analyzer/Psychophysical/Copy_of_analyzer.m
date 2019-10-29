close all;
clear all;

filePath = '/Users/fradeshieh/Desktop/matlab/a3/test';
file = dir(fullfile(filePath,'retinalStabilization_*.mat'));
flipNameStr = 'flip';
pBias = nan(4,length(file));
pThreshold = nan(4,length(file));
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
colorIndex = [0 0 0; 1 0.0 0.0; 0 0.8 0; 0 0 1];

figureNum = 1;
for i = 1:length(file)
    fileName = file(i).name;
    data = load(fullfile(filePath,fileName));
    popChoice = cat(1,popChoice,data.choice);
    trialLength = length(data.Conditions);
    flipDegree = contains(fileName,flipNameStr);
    num = regexp(fileName, '(\d+)','tokens');
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
    
    
    
    for triali = 1:trialLength
        
        if isempty(data.Conditions{triali})
            continue
        else
            motionType = data.Conditions{triali}(end);
        end
        
        if motionType == 1
            headDeg{1}(triali,:) = [data.Conditions{triali}(2),triali];
            headSpe{1}(triali,:) = [data.Conditions{triali}(3),triali];
        elseif motionType == 2
            headDeg{2}(triali,:) = [data.Conditions{triali}(2),triali];
            headSpe{2}(triali,:) = [data.Conditions{triali}(3),triali];
            pursuitDir{2}(triali,:) = [data.Conditions{triali}(4),triali];
            fixSpe{2}(triali,:) = [data.Conditions{triali}(5),triali];
        elseif motionType == 3
            headDeg{3}(triali,:) = [data.Conditions{triali}(2),triali];
            headSpe{3}(triali,:) = [data.Conditions{triali}(3),triali];
            pursuitDir{3}(triali,:) = [data.Conditions{triali}(4),triali];
            rotationSpe{3}(triali,:) = [data.Conditions{triali}(5),triali];
        elseif motionType == 4
            headDeg{4}(triali,:) = [data.Conditions{triali}(2),triali];
            headSpe{4}(triali,:) = [data.Conditions{triali}(3),triali];
            pursuitDir{4}(triali,:) = [data.Conditions{triali}(4),triali];
            fixSpe{4}(triali,:) = [data.Conditions{triali}(5),triali];
        end
    end
    
    uniqueDeg = cell(1,4);
    for typei = 1:4
        del = isnan(headDeg{typei}(:,1));
        headDeg{typei}(del,:) = [];
        
        if flipDegree
            headDeg{typei}(:,1) = -headDeg{typei}(:,1);
        end
        
        del = isnan(headSpe{typei}(:,1));
        headSpe{typei}(del,:) = [];
        if typei == 2 || typei == 4
            del = isnan(pursuitDir{typei}(:,1));
            pursuitDir{typei}(del,:) = [];
            del = isnan(fixSpe{typei}(:,1));
            fixSpe{typei}(del,:) = [];
            
            popPursuitDir{typei} = cat(1,popPursuitDir{typei},pursuitDir{typei});
            popFixSpe{typei} = cat(1,popFixSpe{typei},fixSpe{typei});
        elseif typei ==3
            del = isnan(pursuitDir{typei}(:,1));
            pursuitDir{typei}(del,:) = [];
            del = isnan(rotationSpe{typei}(:,1));
            rotationSpe{typei}(del,:) = [];
            
            popPursuitDir{typei} = cat(1,popPursuitDir{typei},pursuitDir{typei});
            popRotationSpe{typei} = cat(1,popRotationSpe{typei},rotationSpe{typei});
        end
        uniqueDeg{typei} = unique(headDeg{typei}(:,1));
        
        popHeadDeg{typei} = cat(1,popHeadDeg{typei},headDeg{typei});
        popHeadSpe{typei} = cat(1,popHeadSpe{typei},headSpe{typei});
        pChoice{typei} = cat(1,pChoice{typei},data.choice(headDeg{typei}(:,2)));
    end
    
    % [~,hSubplot1] = tight_subplot(1,4,[0.02 0.02]);
    if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
                plot([90,90],[0,1],'-.k');
    for typei=1:4
        choice = cell(length(uniqueDeg{typei}),1);
        pRight = cell(length(uniqueDeg{typei}),1);
        choiceTime = cell(length(uniqueDeg{typei}),1);
        
        if typei > 1
            uniqueDir = unique(pursuitDir{typei}(:,1));
            dirLength = length(uniqueDir);
            for directioni = 1:dirLength
                trialIndexI = find(pursuitDir{typei}(:,1) == uniqueDir(directioni));
                dirIndex{directioni} = pursuitDir{typei}(trialIndexI,:);
            end
        end
        
        if typei == 1
            for k = 1:length(uniqueDeg{typei})
                index = ismember(headDeg{typei}(:,1),uniqueDeg{typei}(k));
                trialIndex = headDeg{typei}(index,2);
                choice{k} = data.choice(trialIndex); %%
                pRight{k} = sum(choice{k} == 2)./ length(choice{k});
                choiceTime{k} = length(choice{k});
            end
            fitData = [uniqueDeg{typei},cell2mat(pRight),cell2mat(choiceTime)];
            
            [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
            pBias(typei,i) = bias;
            pThreshold(typei,i) = threshold;
            xi = min(uniqueDeg{typei}):0.1:max(uniqueDeg{typei});
            y_fit = cum_gaussfit([bias,threshold],xi);
            
            if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
            plot([90,90],[0,1],'-.k');
            hold on
            plot(uniqueDeg{typei},cell2mat(pRight),'*k');
            plot(xi,y_fit,'-k');
            set(gca, 'xlim',[-15,15],'ylim',[0 1])
            xlabel('Heading degree');
            ylabel('Proportion of "right" choice');
            title(['Condition ' num2str(typei) ' for participant ' subNum]);
            %     hleg1=legend('choice','mean & standard error','linear result');
            %     set(hleg1,'Location','EastOutside')
          
           
        elseif typei>1
            for directioni = 1:dirLength
                for k = 1:length(uniqueDeg{typei})
                    index = find(ismember(headDeg{typei}(:,1),uniqueDeg{typei}(k)));
                    trialIndex = headDeg{typei}(index,2);
                    trialIndex = trialIndex(ismember(trialIndex,dirIndex{directioni}(:,2)));
                    choice{k} = data.choice(trialIndex); %%
                    pRight{k} = sum(choice{k} == 2)./ length(choice{k});
                    choiceTime{k} = length(choice{k});
                    
                end
                fitData = [uniqueDeg{typei},cell2mat(pRight),cell2mat(choiceTime)];
                if min(cell2mat(pRight)) == 1
                    fitData(:,2) = fitData(:,2) - 0.01;
                elseif max(cell2mat(pRight)) == 0
                    fitData(:,2) = fitData(:,2) + 0.01;
                end
                [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
                pBias(typei,i) = bias;
                pThreshold(typei,i) = threshold;
                xi = min(uniqueDeg{typei}):0.1:max(uniqueDeg{typei});
                y_fit = cum_gaussfit([bias,threshold],xi);
                
                
                hold on
                plot(uniqueDeg{typei},cell2mat(pRight),'*','color',colorIndex(typei,:));
                plot(xi,y_fit,'-','color',colorIndex(typei,:));
                set(gca, 'xlim',[-15,15],'ylim',[0 1])
                xlabel('Heading degree');
                ylabel('Proportion of "right" choice');
                title(['Condition ' num2str(typei) ' for participant ' subNum]);
                %     hleg1=legend('choice','mean & standard error','linear result');
                %     set(hleg1,'Location','EastOutside')
                 
            end
        end
        figureNum = figureNum+1;
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

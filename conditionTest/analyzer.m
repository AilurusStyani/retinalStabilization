 close all;
clear all;

filePath = fullfile(pwd,'data');
flipNameStr = 'flip';
matFile = dir(fullfile(filePath,'*.mat'));
pBias = cell(1,2);
pThreshold = cell(1,2);
legendText = cell(2,2);

popChoice = [];
pI = 0;

figureNum = 1;
uniqueDeg = cell(1,2);
popHeadDeg = cell(1,2);
popPurDir = cell(1,2);
pChoice = cell(1,2);

colorIndex = {[1 0.1 0.1],[0.8 0.5 0.3];[0.1 0.1 1],[0.3 0.5 0.8]}; % {typeI,dir}

for i = 1:length(matFile)
    
    data = load(fullfile(filePath,matFile(i).name));
    popChoice = cat(1,popChoice,data.choice);
    trialNum = size(data.conditionIndex,1);
    flipDegree = contains(matFile(i).name,flipNameStr);
    num = regexp(matFile(i).name, '(\d+)','tokens');
    subNum = cell2mat(num{1});
    
    % rotation:
    headDeg{1} = nan(trialNum,2);
    purDir{1} = nan(trialNum,2);
    choice{1} = nan(trialNum,2);
    
    % frustum shift
    headDeg{2} = nan(trialNum,2);
    purDir{2} = nan(trialNum,2);
    choice{2} = nan(trialNum,2);
    
    % classfy the degree and pursuit direction based on the rotation type
    for triali = 1:trialNum
        motionType = data.trialIndex(triali,2);
        headDeg{motionType}(triali,:) = [data.conditionIndex(triali,1),data.trialIndex(triali,1)];
        purDir{motionType}(triali,:) = [data.conditionIndex(triali,4),data.trialIndex(triali,1)];
        choice{motionType}(triali,:) = [data.choice(triali),data.trialIndex(triali,1)];
    end
    
    if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white'); hold on; plot([90,90],[0,1],'-.k');
    
    for typeI =1:length(headDeg)
        % delete the nan vaule
        delDeg = find(isnan(headDeg{typeI}(:,1)));
        headDeg{typeI}(delDeg,:) = [];
        delDir = find(isnan(purDir{typeI}(:,1)));
        purDir{typeI}(delDir,:) = [];
        delChoice = find(isnan(choice{typeI}(:,1)));
        choice{typeI}(delChoice,:) = [];
        
        if flipDegree
            headDeg{typeI}(:,1) = -headDeg{typeI}(:,1);
        end
        
        uniqueDeg{typeI} = unique(headDeg{typeI}(:,1));
        
        popHeadDeg{typeI} = cat(1,popHeadDeg{typeI},headDeg{typeI});
        popPurDir{typeI} = cat(1,popPurDir{typeI},purDir{typeI});
        pChoice{typeI} = cat(1,pChoice{typeI},choice{typeI});
        
        pRight = [];
        choiceF = [];
        
        uniqueDir = unique(purDir{typeI}(:,1));
        dirLength = length(uniqueDir);
        dirIndex = cell(1,dirLength);
        for dir = 1:dirLength
            trialIndexI = find(purDir{typeI}(:,1) == uniqueDir(dir));
            dirIndex{dir} = purDir{typeI}(trialIndexI,:);
        end
        
        %         if typeI == 1
        %             for k = 1:length(uniqueDeg{typeI})
        %                 index = ismember(headDeg{typeI}(:,1),uniqueDeg{typeI}(k));
        %                 trialIndex = headDeg{typeI}(index,2);
        %                 pRight{k} = sum(choice{k} == 2)./ length(choice{k});
        %                 choiceF{k} = length(choice{k});
        %             end
        %             fitData = [uniqueDeg{typeI},cell2mat(pRight),cell2mat(choiceF)];
        %
        %             [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
        %             pBias(triali,i) = bias;
        %             pThreshold(triali,i) = threshold;
        %             xi = min(uniqueDeg{typeI}):0.1:max(uniqueDeg{typeI});
        %             y_fit = cum_gaussfit([bias,threshold],xi);
        %
        %             if ishandle(figureNum); end; figure(figureNum); set(gcf,'color','white');
        %             plot([90,90],[0,1],'-.k');
        %             hold on
        %             plot(uniqueDeg{typeI},cell2mat(pRight),'*');
        %             plot(xi,y_fit,'-');
        %             set(gca, 'xlim',[-15,15],'ylim',[0 1])
        %             xlabel('Heading degree');
        %             ylabel('Proportion of "right" choice');
        %             title(['Condition ' num2str(typeI) ' for participant ' subNum]);
        %             %     hleg1=legend('choice','mean & standard error','linear result');
        %             %     set(hleg1,'Location','EastOutside')
        %             text(5,0.8,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color','b')
        %             text(5,0.7,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color','b')
        %             figureNum = figureNum+1;
        %         elseif typeI>1
        for dir = 1:dirLength
            for k = 1:length(uniqueDeg{typeI})
                index = find(ismember(headDeg{typeI}(:,1),uniqueDeg{typeI}(k)));
                trialIndex = headDeg{typeI}(index,2);
                trialIndex = trialIndex(ismember(trialIndex,dirIndex{dir}(:,2)));
                choicei{k} = choice{typeI}(ismember(choice{typeI}(:,2),trialIndex),1); %%
                pRight{k} = sum(choicei{k} == 2)./ length(choicei{k});
                choiceF{k} = length(choicei{k});
            end
            fitData = [uniqueDeg{typeI},cell2mat(pRight'),cell2mat(choiceF')];
            if min(cell2mat(pRight)) == 1
                fitData(:,2) = fitData(:,2) - 0.01;
            elseif max(cell2mat(pRight)) == 0
                fitData(:,2) = fitData(:,2) + 0.01;
            end
            [bias,threshold] = cum_gaussfit_max1(fitData(1:end,:));
            pBias{typeI}(dir,i) = bias;
            pThreshold{typeI}(dir,i) = threshold;
            xi = min(uniqueDeg{typeI}):0.1:max(uniqueDeg{typeI});
            y_fit = cum_gaussfit([bias,threshold],xi);
            
            plot(uniqueDeg{typeI},cell2mat(pRight),'*','color',colorIndex{dir,typeI});
            h(typeI,dir) = plot(xi,y_fit,'-','color',colorIndex{dir,typeI});
            set(gca, 'xlim',[-15,15],'ylim',[0 1])
            xlabel('Heading degree');
            ylabel('Proportion of "right" choice');
            legendText{dir,typeI} = ['Condition ' num2str(typeI) ' direction ' num2str(dir)];
            %     hleg1=legend('choice','mean & standard error','linear result');
            %     set(hleg1,'Location','EastOutside')
            text(6*typeI,0.6-0.2*dir,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color',colorIndex{dir,typeI})
            text(6*typeI,0.5-0.2*dir,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color',colorIndex{dir,typeI})
        end
        %         end
    end
    legend(h(:),cell2mat(legendText(:)));
    title([ 'Participant''s Number: ' subNum]);
    figureNum = figureNum+1;
    
end

if ishandle(figureNum); end; figure(figureNum); clf; set(gcf,'color','white'); hold on; plot([90,90],[0,1],'-.k');
% hPopulation = tight_subplot(2,2,[0.1 0.1]);
% set(hPopulation,'xtick',[]);
% set(hPopulation,'ytick',[]);

% draw popular result
for typeI =1:length(popHeadDeg)
    pPRight = [];
    pChoiceF = [];
    
    popUniqueDeg{typeI} = unique(popHeadDeg{typeI}(:,1));
    pUniqueDir = unique(popPurDir{typeI}(:,1));
    dirLength = length(pUniqueDir);
    dirIndex = cell(1,dirLength);
    for dir = 1:dirLength
        trialIndexI = find(popPurDir{typeI}(:,1) == pUniqueDir(dir));
        dirIndex{dir} = popPurDir{typeI}(trialIndexI,:);
    end
    for dir = 1:dirLength
        for k = 1:length(popUniqueDeg{typeI})
            index = find(ismember(popHeadDeg{typeI}(:,1),popUniqueDeg{typeI}(k)));
            trialIndex = popHeadDeg{typeI}(index,2);
            trialIndex = trialIndex(ismember(trialIndex,dirIndex{dir}(:,2)));
            pChoicei{k} = pChoice{typeI}(ismember(pChoice{typeI}(:,2),trialIndex),1); %%
            pPRight{k} = sum(pChoicei{k} == 2)./ length(pChoicei{k});
            pChoiceF{k} = length(pChoicei{k});
        end
        pFitData = [popUniqueDeg{typeI},cell2mat(pPRight'),cell2mat(pChoiceF')];
        if min(cell2mat(pPRight)) == 1
            pFitData(:,2) = pFitData(:,2) - 0.01;
        elseif max(cell2mat(pPRight)) == 0
            pFitData(:,2) = pFitData(:,2) + 0.01;
        end
        [bias,threshold] = cum_gaussfit_max1(pFitData(1:end,:));
        xi = min(popUniqueDeg{typeI}):0.1:max(popUniqueDeg{typeI});
        y_fit = cum_gaussfit([bias,threshold],xi);
        
        plot(popUniqueDeg{typeI},cell2mat(pPRight),'*','color',colorIndex{dir,typeI});
        h(typeI,dir) = plot(xi,y_fit,'-','color',colorIndex{dir,typeI});
        set(gca, 'xlim',[-15,15],'ylim',[0 1])
        xlabel('Heading degree');
        ylabel('Proportion of "right" choice');
        legendText{dir,typeI} = ['Condition ' num2str(typeI) ' direction ' num2str(dir)];
        %     hleg1=legend('choice','mean & standard error','linear result');
        %     set(hleg1,'Location','EastOutside')
        text(6*typeI,0.6-0.2*dir,sprintf('\\it\\mu_{psy} = \\rm%6.3g\\circ',bias),'color',colorIndex{dir,typeI})
        text(6*typeI,0.5-0.2*dir,sprintf('\\it\\sigma_{psy} = \\rm%6.3g\\circ', threshold),'color',colorIndex{dir,typeI})
    end
    %         end
end
legend(h(:),cell2mat(legendText(:)));
title('Population result');
[h_Bias,p_Bias] = ttest2(pBias{1},pBias{2})
[h_Threshold,p_Threshold] = ttest2(pThreshold{1},pThreshold{2})

% originally modified by GY from GCD (cum_gaussfit_max.m)
function [alpha, beta] = cum_gaussfit_max1(data_cum)
% Gaussian fits a accumulative Gaussian function to data using maximum likelihood
%	maximization under binomial assumptions.  It uses Gaussian_Fun for
%	error calculation.  Data must be in 3 columns: x, %-correct, Nobs
%	usage: [alpha, beta] = cum_gaussfit_max1(data)
%		alpha and beta are the bias and standard deviation parameters.
global Data_cum;
Data_cum = data_cum;

% generate guess
q = ones(2,1);

% get a starting threshold estimate by testing a bunch of values
% bias range from [-100,-10,-1,0,1,10,100]
% threshold ranges from [0.1,1,10,100]
bias_e = [-100,-10,-1,0,1,10,100];
threshold_e = [0.1,1,10,100];
errors=[];
for i=1:length(bias_e)
    for j=1:length(threshold_e)
        q(1,1) = bias_e(i);
        q(2,1) = threshold_e(j);
        errors(i,j) = cum_gaussfit_max(q);
    end
end
[min_indx1,min_indx2] = find(errors==min(min(errors)));
q(1,1) = bias_e(min_indx1);
q(2,1) = threshold_e(min_indx2);

OPTIONS = optimset('MaxIter', 5000);
quick = fminsearch('cum_gaussfit_max',q);
% quick
alpha = quick(1,1);
beta = quick(2,1);
end

function p=cum_gaussfit(beta,x)
MU=beta(1);
SIGMA=beta(2);
p = normcdf(x,MU,SIGMA);
end
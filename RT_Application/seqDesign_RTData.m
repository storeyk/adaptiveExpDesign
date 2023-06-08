
%This function runs the sequential design algorithm described in Cho,
%Lewis, Storey, Zittle, "An adaptive information-theoretic...", 2023,
%using the synthetic RT data generated from the hybrid CA model.

%   Inputs:
%       -responderType: uses CA data - enter 'High','Low', or 'Med'
%               -'High': tau = 15, alpha = .14, beta = .093
%               -'Med': tau = 22, alpha = .05, beta = .033
%               -'Low': tau = 22, alpha = .014, beta = .014
%       -modelType: choices are 'EXP+DVR','LOG+DVR','LOG+CCR', as described in
%               Mohsin, Enderling et. al, 2002
%
%   Note: This script requires use of Marko Laine's MCMC toolbox code
%   (found at https://mjlaine.github.io/mcmcstat/#orgcdeadeb)
%
% Ref: TBD
%
% Author: Allison Lewis <lewisall@lafayette.edu>
% Last revision: 06-07-2023

function seqDesign_RTData(responderType, getmodelType)


global currentPts modelType data growthParams

modelType = getmodelType; 

% First, load data
filename = ['Data/CAData_' responderType 'Responder.mat'];
load(filename)
fulldata.xdata = data.xdata(1:55);
fulldata.ydata = data.ydata(1:55,1); %use only tumor volume during treatment period
clear data


path = ['Figures/' responderType 'Responder/' modelType '/'];
mkdir( path )


%Settings
no_smps = 10000; %Metropolis chain lengths

%Set up vectors for saved quantities;
beta_vec = []; %store parameter estimates at each iteration
err_vec = []; %store MSE at each iteration
widthCI = []; %store width of credible interval at each iteration
areaCI = []; %store area of credible interval at each iteration
penCoeff_vec = [];

%Initial starting data - days 8 & 15 were used to get intrinsic growth
%params, then add day 19 to start the RT parameter estimation
currentPts = [fulldata.xdata(9) fulldata.ydata(9);
    fulldata.xdata(16) fulldata.ydata(16);
    fulldata.xdata(20) fulldata.ydata(20)];


%Possible experimental designs - discrete times in future, starting
%w/day 20
expDesigns = fulldata.xdata((21):end);


%% Set intrinsic growth params using two pre-treatment data points
% (calibrations done in separate script)

switch modelType
    case 'EXP+DVR'
        switch responderType
            case {'High'}
                growthParams = .1914;
            case {'Med'}
                growthParams = .2103;
            case {'Low'}
                growthParams = .2081;
        end
    case {'LOG+DVR','LOG+CCR'}
        switch responderType
            case {'High'}
                growthParams = [.3178 .3984];
            case {'Med'}
                growthParams = [.3952 .4782];
            case {'Low'}
                growthParams = [.4148 .4527];
        end
end



%% Now start the sequential design procedure, estimating beta only

nIter = 0; %keep track of how many calibration iterations have been performed
ptsLeft = 1; %binary tracker; 1 means there are still points left to choose

switch modelType
    case {'EXP+DVR','LOG+DVR'}
        beta_val = .1; %initial guess
    case {'LOG+CCR'}
        beta_val = 0.001; %initial guess
end


while ptsLeft == 1
    
    nIter = nIter+1;
    
    params1 = {
        {'beta',beta_val,0,1}
        };
    
    data.xdata = currentPts(:,1);
    data.ydata = currentPts(:,2);
    
    model.ssfun = @ssq_tumorVolpostRT;
    options.updatesigma = 1;
    
    options.nsimu = no_smps;
    [results,chain,s2chain] = mcmcrun(model,data,params1,options);
    
    options.nsimu = no_smps;
    [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);
    
    
    % Find the optimal parameters
    ind = find(ss2chain == min(ss2chain));  ind = ind(1);
    beta_val = chain(ind,:); %This is our fitted parameter set
    params = [growthParams beta_val];
    
    %Save off parameter estimates and chains
    beta_vec = [beta_vec; beta_val];
    fitChain{nIter}.paramschain = chain;
    fitChain{nIter}.s2chain = s2chain;
    fitChain{nIter}.ss2chain = ss2chain;
    
    
    % Generate current model trajectory, credible intervals, and plot
    [timeFit, volFit] = tumorModel_postRT(modelType,params);
    lowfi = interp1(timeFit,volFit, fulldata.xdata);
    
    ln = fulldata.xdata';
    modelfun = @(ln,params) tumorfun(ln,params);
    
    figure(1)
    pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
    pred.obslims = [];
    mcmcpredplot(pred)
    hold on
    h=gca;
    plot(fulldata.xdata,fulldata.ydata,'ok','MarkerSize',6) %All possible scans
    plot(fulldata.xdata,lowfi,'-k','Linewidth',2) %This is your optimal model fit
    plot(currentPts(:,1),currentPts(:,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
    xlabel('Time (days)','FontSize',18);
    ylabel('Volume','FontSize',18);
    titleName = [responderType 'Responder: Iteration ' num2str(nIter)];
    title(string(titleName),'Interpreter','none','FontSize',14)
    axis([min(fulldata.xdata)-1 max(fulldata.xdata)+1 0 max(fulldata.ydata)+.2])
    set(h, 'FontSize',18)
    hold off
    filename = [path 'ModelFit_Iteration' num2str(nIter) '.jpg'];
    saveas(gcf,filename)
    filename = [path 'ModelFit_Iteration' num2str(nIter) '.fig'];
    saveas(gcf,filename)
    close(gcf)
    
    
    figure(2)
    mcmcplot(chain,[],{'\beta'},'chainpanel')
    filename = [path 'Chain_Iteration' num2str(nIter) '.jpg'];
    saveas(gcf,filename)
    filename = [path 'Chain_Iteration' num2str(nIter) '.fig'];
    saveas(gcf,filename)
    close(gcf)
    
    figure(3)
    mcmcplot(chain,[],{'\beta'},'denspanel')
    filename = [path 'Density_Iteration' num2str(nIter) '.jpg'];
    saveas(gcf,filename)
    filename = [path 'Density_Iteration' num2str(nIter) '.fig'];
    saveas(gcf,filename)
    close(gcf)
    
    
    % Calculate MSE to measure error
    mse = sum((lowfi'-fulldata.ydata).^2)/numel(fulldata.ydata);
    err_vec = [err_vec; mse];
    
    % Calculate uncertainty metric
    for n = 1:length(ln)
        widthCI(nIter,n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
    end
    areaCI = [areaCI; sum(widthCI(nIter,1:end))];
    
    
    
    % Determine which point should be selected next
    if currentPts(end,1)==fulldata.xdata(end)
        ptsLeft = 0; %this is the final iteration
        
    else
        
        knnSet = 10:10:no_smps; %thin out samples for kNN analysis
        newChain = chain(knnSet,:);
        
        % Calculate predicted values at each experimental design with
        % each parameter set from newChain
        for ii = 1:size(newChain,1)
            
            [tsol, ysol] = tumorModel_postRT(modelType,[growthParams newChain(ii,:)]);
            
            for jj = 1:length(expDesigns)
                lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
            end
        end
        
        %Calculate MI for each remaining design        

        % Standardize chain and low-fi predictions to N(0,1)
        normChain = (newChain-mean(newChain))./std(newChain);
        normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);
        
        for jj = 1:length(expDesigns)
            [I1,~] = KraskovMI(normChain,normlowfiOut(:,jj),6); 
            miVals(jj) = I1;
        end
        maxMI = max(miVals);
        relMI = miVals/maxMI;
        lbMI = min(miVals);
        ubMI = max(miVals);
        if lbMI==ubMI
            unifMI = ones(1,length(miVals));
        else
            unifMI = (miVals-lbMI)/(ubMI-lbMI);
        end
        
        
        %Compute penalty coefficient
        pen = abs(lowfi(end)-currentPts(end,2))/(lowfi(end)+currentPts(end,2));
        penCoeff_vec = [penCoeff_vec pen];
        
        %Penalize MI values to get score function
        for p = 1:length(relMI)
            score(p) = unifMI(p) - abs(pen)*sum(unifMI(1:(p-1)))/sum(unifMI);
        end

        figure(1000)
        plot(expDesigns,score,'ob','MarkerFaceColor','b')
        hold on
        plot(expDesigns,unifMI,'or')
        filename = [path 'ScoreFxn_Iteration' num2str(nIter) '.jpg'];
        saveas(gcf,filename)
        close(gcf)

        
        %Choose optimal design
        point = expDesigns(1,find(score == max(score)))
        
        if length(point)>1
            fprintf('Multiple points have same score function. Last point in tie list chosen.')
            point = point(end);
        end
        
        
        %Add point to current list of points for next round of calibration
        idx = find(fulldata.xdata==point);
        currentPts = [currentPts; point fulldata.ydata(idx)]
        
        %Remove chosen point and all skipped points from experimental design list
        expDesigns(find(expDesigns<=point)) = [];
        
        
        clear score relMI miVals
    end
end



point_list = currentPts; %Final list of points in order of selection


% Save results for later comparison
file = [path 'Results.mat'];
save(file, 'point_list','beta_vec','growthParams','err_vec','penCoeff_vec','areaCI','fitChain')


end



%% Helper functions


% SSQ function for calibration of RT parameter
function SSrt = ssq_tumorVolpostRT(params,data)

global currentPts modelType growthParams

[time,vol] = tumorModel_postRT(modelType,[growthParams params]);

tumVol = interp1(time,vol, currentPts(:,1));

SSrt = sum((tumVol - currentPts(:,2)).^2);

end

% Function for credible interval plotting
function v=tumorfun(timef,params)

global modelType growthParams


[tsol, ysol] = tumorModel_postRT(modelType,[growthParams params]);

v = interp1(tsol,ysol,timef);

end




%This function runs the sequential design algorithm described in Cho,
%Lewis, Storey, Zittle, "An adaptive information-theoretic...", 2023,
%using the IAS therapy data from Bruchovsky et. al (2009).
%
%   Note: This script requires use of Marko Laine's MCMC toolbox code
%   (found at https://mjlaine.github.io/mcmcstat/#orgcdeadeb)
%
% Ref: TBD
%
% Author: Allison Lewis <lewisall@lafayette.edu>
% Last revision: 06-07-2023

function seqDesign_IASData_PSAonly

global currentPts data fixedVals
 

% First, load data
filename = ['Data/Patient39Data.mat'];
load(filename)
dataIdx = ~isnan(patientData.psa);
numPts = sum(dataIdx);
fulldata.xdata = patientData.day(dataIdx); %currently set for all cycles
fulldata.ydata = patientData.psa(dataIdx); %use psa as y data
fulldata.xdata = fulldata.xdata(1:24); 
fulldata.ydata = fulldata.ydata(1:24); %end at completion of second round of treatment (1.5 cycles) 
Ainit=patientData.testosterone(1);
Pinit=patientData.psa(1);


path = ['Figures/Patient39/PSAonly/1.5Cycles/'];
mkdir( path )
%Settings
no_smps = 10000; %Metropolis chain lengths

%Set up vectors for saved quantities;
q2_vec = [];
gamma1_vec = []; 
sigma2_vec = [];
A0_vec = [];

err_vec = []; %store MSE at each iteration
widthCI = []; %store width of credible interval at each iteration
areaCI = []; %store area of credible interval at each iteration
penCoeff_vec = [];

%Initial starting data - using points 1&2 (assuming other parameters already
%fixed in advance
currentPts = [fulldata.xdata(1:2) fulldata.ydata(1:2)];


%Possible experimental designs - discrete times in future
expDesigns = fulldata.xdata(3:end);

%% Fix non-treatment parameters (estimates found using fmincon)
load('paramsFixedValues.mat')


%% Now start the sequential design procedure, estimating [q2, gamma1, A0, sigma2]

nIter = 0; %keep track of how many calibration iterations have been performed
ptsLeft = 1; %binary tracker; 1 means there are still points left to choose

initGuess = fixedVals;
[lb,ub,~,~] = getParamInfo;

while ptsLeft == 1
    
    nIter = nIter+1;
    
    params1 = {
        {'q2',initGuess(3),lb(3),ub(3)}
        {'gamma_1',initGuess(7),lb(7),ub(7)}
        {'A0',initGuess(9),lb(9),ub(9)}
        {'sigma2',initGuess(14),lb(14),ub(14)}
        };
    
    data.xdata = currentPts(:,1);
    data.ydata = currentPts(:,2);
    
    model.ssfun = @(params,data)ssq_PCmodel(params,fulldata,dayswitch,Ainit,Pinit);
    options.updatesigma = 1;
    
    options.nsimu = no_smps/5;
    [results,chain,s2chain] = mcmcrun(model,data,params1,options);
    
    options.nsimu = no_smps;
    [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);
    
    
    % Find the optimal parameters
    estParams = mean(chain); 
    estParams = results.mean;
    params = [fixedVals(1:2) estParams(1) fixedVals(4:6) estParams(2) fixedVals(8) estParams(3) fixedVals(10:13) estParams(4) fixedVals(15:18)];

    %Save off parameter estimates and chains
    q2_vec = [q2_vec; estParams(1)];
    gamma1_vec = [gamma1_vec; estParams(2)];
    A0_vec = [A0_vec; estParams(3)];
    sigma2_vec = [sigma2_vec; estParams(4)];
    fitChain{nIter}.paramschain = chain;
    fitChain{nIter}.s2chain = s2chain;
    fitChain{nIter}.ss2chain = ss2chain;
    
    
    % Generate current model trajectory, credible intervals, and plot
    [~,lowfi] = evaluatePCmodel(params,fulldata,dayswitch,Ainit,Pinit);

    ln = (fulldata.xdata(1):10:fulldata.xdata(end))';
    modelfun = @(ln,params) predPCmodel(ln,params,fulldata,dayswitch,Ainit,Pinit);
    
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
    ylabel('PSA','FontSize',18);
    titleName = ['Patient 39: Iteration ' num2str(nIter)];
    title(string(titleName),'Interpreter','none','FontSize',14)
    set(h, 'FontSize',18)
    hold off
    filename = [path 'ModelFit_Iteration' num2str(nIter) '.jpg'];
    saveas(gcf,filename)
     filename = [path 'ModelFit_Iteration' num2str(nIter) '.fig'];
     saveas(gcf,filename)
    close(gcf)
    
    
    figure(2)
    mcmcplot(chain,[],{'q_2','\gamma_1','A_0','\sigma_2'},'chainpanel')
    filename = [path 'Chain_Iteration' num2str(nIter) '.jpg'];
    saveas(gcf,filename)
    filename = [path 'Chain_Iteration' num2str(nIter) '.fig'];
    saveas(gcf,filename)
    close(gcf)
    
    figure(3)
    mcmcplot(chain,[],{'q_2','\gamma_1','A_0','\sigma_2'},'denspanel')
    filename = [path 'Density_Iteration' num2str(nIter) '.jpg'];
    saveas(gcf,filename)
    filename = [path 'Density_Iteration' num2str(nIter) '.fig'];
    saveas(gcf,filename)
    close(gcf)
    
    
    % Calculate MSE to measure error
    mse = sum((lowfi-fulldata.ydata).^2)/numPts;
    err_vec = [err_vec; mse];
    
    %Calculate uncertainty metrics
    for n = 1:length(ln)
        widthCI(nIter,n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
    end
    %multiply each CIwidth by 10 because we're using time steps of 10 days
    areaCI = [areaCI; sum(10*widthCI(nIter,1:end))];

    
    % Determine which point should be selected next
    if currentPts(end,1)==fulldata.xdata(end)
        ptsLeft = 0; %this is the final iteration
        
    else
        
        knnSet = 10:10:no_smps; %thin out samples for kNN analysis
        newChain = chain(knnSet,:);
        
        % Calculate predicted values at each experimental design with
        % each parameter set from newChain
        for ii = 1:size(newChain,1)
            
            [tsol, ysol] = evaluatePCmodel([fixedVals(1:2) newChain(ii,1) fixedVals(4:6) newChain(ii,2) fixedVals(8) newChain(ii,3) fixedVals(10:13) newChain(ii,4) fixedVals(15:18)],fulldata,dayswitch,Ainit,Pinit);
            
            for jj = 1:length(expDesigns)
                lowfiOut(ii,jj) = interp1(tsol,ysol,expDesigns(jj));
            end
        end
        
        %Calculate MI for each remaining design        
     
        %Standardize to N(0,1)
        normChain = (newChain-mean(newChain))./std(newChain);
        normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);
        
        for jj = 1:length(expDesigns)
            [I1,~] = KraskovMI(normChain,normlowfiOut(:,jj),6); 
            miVals(jj) = I1;
        end
        maxMI = max(miVals);
        lbMI = min(miVals);
        ubMI = max(miVals);
        if lbMI==ubMI
            unifMI = ones(1,length(miVals));
        else
            unifMI = (miVals-lbMI)/(ubMI-lbMI);
        end
        
        %Calculate score function by penalizing for points skipped
        pen = abs(lowfi(end)-currentPts(end,2))/(lowfi(end)+currentPts(end,2));
        penCoeff_vec = [penCoeff_vec pen];
        
        %penalize MI by k*absolute rcg*penalty for skipped points
        for p = 1:length(unifMI)
            score(p) = unifMI(p) - abs(pen)*sum(unifMI(1:(p-1)))/sum(unifMI);
        end

        figure(1000)
        plot(expDesigns,score,'ob','MarkerFaceColor','b')
        hold on
        plot(expDesigns,unifMI,'or')
        filename = [path 'ScoreFxn_Iteration' num2str(nIter) '.fig'];
        saveas(gcf,filename)
        filename = [path 'ScoreFxn_Iteration' num2str(nIter) '.jpg'];
        saveas(gcf,filename)
        close(gcf)

        
        %Choose optimal design
        point = expDesigns(find(score == max(score)))
        
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
save(file, 'point_list','q2_vec','gamma1_vec','A0_vec','sigma2_vec','fixedVals','err_vec','penCoeff_vec','areaCI','fitChain')


end



%% Helper functions


% SSQ function for calibration of RT parameter
function SSrt = ssq_PCmodel(params,data,dayswitch,Ainit,Pinit)

global currentPts fixedVals

paramSet = [fixedVals(1:2) params(1) fixedVals(4:6) params(2) fixedVals(8) params(3) fixedVals(10:13) params(4) fixedVals(15:18)];

[time,psamodel] = evaluatePCmodel(paramSet,data,dayswitch,Ainit,Pinit);

psaoutput = interp1(time,psamodel, currentPts(:,1));

SSrt = sum((psaoutput - currentPts(:,2)).^2);

end

% Function for credible interval plotting
function v=predPCmodel(timef,params,data,dayswitch,Ainit,Pinit)

global fixedVals

paramSet = [fixedVals(1:2) params(1) fixedVals(4:6) params(2) fixedVals(8) params(3) fixedVals(10:13) params(4) fixedVals(15:18)];


[tsol, ysol] = evaluatePCmodel(paramSet,data,dayswitch,Ainit,Pinit);

v = interp1(tsol,ysol,timef);

end



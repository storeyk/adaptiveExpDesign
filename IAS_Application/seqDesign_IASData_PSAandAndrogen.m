
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

function seqDesign_IASData_PSAandAndrogen


global currentPts data fixedVals
 

% First, load data
filename = ['Data/Patient39Data.mat'];
load(filename)
dataIdx = ~isnan(patientData.psa);
numPts = sum(dataIdx);
fulldata.xdata = patientData.day(dataIdx); 
fulldata.psadata = patientData.psa(dataIdx); 
fulldata.anddata = patientData.testosterone(dataIdx); 
fulldata.xdata = fulldata.xdata(1:24); 
fulldata.psadata = fulldata.psadata(1:24); %24 - end at completion of second round of treatment (1.5 cycles) 60 - end at completion of a treatment cycle
fulldata.anddata = fulldata.anddata(1:24);
Ainit=patientData.testosterone(1);
Pinit=patientData.psa(1);


path = ['Figures/Patient39/BothData/1.5Cycles/'];
mkdir( path )


%Settings
no_smps = 10000; %Metropolis chain lengths

%Set up vectors for saved quantities;
q2_vec = [];
gamma1_vec = []; %store parameter estimates at each iteration
sigma2_vec = [];
A0_vec = [];

err_vec = []; %store MSE at each iteration
widthCI = []; %store width of credible interval at each iteration
areaCI = []; %store area of credible interval at each iteration
pen_vec = [];
penCoeffPSA_vec = [];
penCoeffAND_vec = [];


%%
%Initial starting data - using first two PSA and first two androgen points
%(attach data type indicator at end: 1 for PSA, 2 for androgen)
currentPts = [fulldata.xdata(1:2) fulldata.psadata(1:2) repmat(1,2,1);
    fulldata.xdata(1:2) fulldata.anddata(1:2) repmat(2,2,1)];

%Possible experimental designs - discrete times in future
expDesigns(1,:) = [fulldata.xdata(3:end)' fulldata.xdata(3:end)'];
expDesigns(2,:) = [ones(1,length(fulldata.xdata(3:end))) 2*ones(1,length(fulldata.xdata(3:end)))];

%% Fix non-treatment parameters (estimates based on fmincon)
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
    
    model.ssfun = @(params,data)ssq_PCmodel_BothData(params,fulldata,dayswitch,Ainit,Pinit);
    options.updatesigma = 1;
    
    options.nsimu = no_smps/5;
    [results,chain,s2chain] = mcmcrun(model,data,params1,options);
    
    options.nsimu = no_smps;
    [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);
    
    
    % Find the optimal parameters
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
    [~,lowfi] = evaluatePCmodel_BothData(params,fulldata,dayswitch,Ainit,Pinit);
    psaIdx = find(currentPts(:,3)==1);
    andIdx = find(currentPts(:,3)==2);
    
    
    ln = (fulldata.xdata(1):10:fulldata.xdata(end))';
    modelfun = @(ln,params) predPCmodel_BothData(ln,params,fulldata,dayswitch,Ainit,Pinit);
    
    figure(1)
    pred = mcmcpred(results,chain,s2chain,ln,modelfun,no_smps);
    pred.obslims = [];
    mcmcpredplot(pred)
    subplot(2,1,1)
    hold on
    h=gca;
    plot(fulldata.xdata,fulldata.psadata,'ok','MarkerSize',6) %All possible scans
    plot(fulldata.xdata,lowfi(:,1),'-k','Linewidth',2) %This is your optimal model fit
    plot(currentPts(psaIdx,1),currentPts(psaIdx,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
    xlabel('Time (days)','FontSize',18);
    ylabel('PSA','FontSize',18);
    title('');
    set(h, 'FontSize',18)
    hold off
    subplot(2,1,2)
    hold on
    h=gca;
    plot(fulldata.xdata,fulldata.anddata,'ok','MarkerSize',6) %All possible scans
    plot(fulldata.xdata,lowfi(:,2),'-k','Linewidth',2) %This is your optimal model fit
    plot(currentPts(andIdx,1),currentPts(andIdx,2),'ok','MarkerFaceColor','k','Linewidth',2,'MarkerSize',8) %Selected scans
    xlabel('Time (days)','FontSize',18);
    ylabel('Androgen','FontSize',18);
    title('');
    set(h, 'FontSize',18)
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
    mse = sum((lowfi(:,1)-fulldata.psadata).^2)/length(fulldata.psadata)+sum((lowfi(:,2)-fulldata.anddata).^2)/length(fulldata.anddata);
    err_vec = [err_vec; mse];
    
    %Calculate uncertainty metrics
    for n = 1:length(ln)
        widthCIpsa(nIter,n) = pred.predlims{1,1}{1}(3,n)-pred.predlims{1,1}{1}(1,n);
        widthCIand(nIter,n) = pred.predlims{1,1}{2}(3,n)-pred.predlims{1,1}{2}(1,n);
    end
    %multiply each CIwidth by 10 because we're using time steps of 10 days
    areaCI = [areaCI; sum(10*widthCIpsa(nIter,1:end))+sum(10*widthCIand(nIter,1:end))];
    
    
    
    % Determine which point should be selected next
    if (length(expDesigns(1,:))==0)
        ptsLeft = 0; %this is the final iteration
        
    else
        
        knnSet = 10:10:no_smps; %thin out samples for kNN analysis
        newChain = chain(knnSet,:);
        
        % Calculate predicted values at each experimental design with
        % each parameter set from newChain
        for ii = 1:size(newChain,1)
            
            [tsol, ysol] = evaluatePCmodel_BothData([fixedVals(1:2) newChain(ii,1) fixedVals(4:6) newChain(ii,2) fixedVals(8) newChain(ii,3) fixedVals(10:13) newChain(ii,4) fixedVals(15:18)],fulldata,dayswitch,Ainit,Pinit);
            
            for jj = 1:length(expDesigns(1,:))
                
                %Pull either tumor volume or necrotic data depending on
                %experimental design
                if expDesigns(2,jj)==1 %interpolate PSA solution
                    lowfiOut(ii,jj) = interp1(tsol, ysol(:,1), expDesigns(1,jj));
                else %interpolate androgen solution
                    lowfiOut(ii,jj) = interp1(tsol, ysol(:,2), expDesigns(1,jj));
                end
            end
        end
        
        %Calculate MI for each remaining design        
     
        % Standardize to N(0,1)
        normChain = (newChain-mean(newChain))./std(newChain);
        normlowfiOut = (lowfiOut - mean(lowfiOut))./std(lowfiOut);
        
        for jj = 1:length(expDesigns(1,:))
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
        
        
        %Penalize using the larger of the two penalty coefficients from multiple metrics
        lastPSA = find(currentPts(:,3)==1,1,'last');
        penCoeffPSA = abs(lowfi(end,1)-currentPts(lastPSA,2))/(lowfi(end,1)+currentPts(lastPSA,2));
        penCoeffPSA_vec = [penCoeffPSA_vec penCoeffPSA];
        lastAND = find(currentPts(:,3)==2,1,'last');
        if isempty(lastAND)==0
            penCoeffAND = abs(lowfi(end,2)-currentPts(lastAND,2))/(lowfi(end,2)+currentPts(lastAND,2));
        else
            penCoeffAND = 0;
        end
        penCoeffAND_vec = [penCoeffAND_vec penCoeffAND];
        pen = max(penCoeffPSA,penCoeffAND); %choose the larger of the two to use in penalty term
        pen_vec = [pen_vec pen];
        
        %penalize MI 
        for p = 1:length(unifMI)
            timePt = expDesigns(1,p);
            skipPts = find(expDesigns(1,:)<timePt);
            score(p) = unifMI(p) - abs(pen)*sum(unifMI(skipPts))/sum(unifMI);
        end

        figure(1000)
        PSAidx = find(expDesigns(2,:)==1);
        ANDidx = find(expDesigns(2,:)==2);
        plot(expDesigns(1,PSAidx),score(PSAidx),'ob','MarkerFaceColor','b')
        hold on
        plot(expDesigns(1,ANDidx),score(ANDidx),'or','MarkerFaceColor','r')
        filename = [path 'ScoreFxn_Iteration' num2str(nIter) '.fig'];
        saveas(gcf,filename)
        filename = [path 'ScoreFxn_Iteration' num2str(nIter) '.jpg'];
        saveas(gcf,filename)
        close(gcf)

     
        %Choose optimal design
        optIdx = find(score==max(score));
        if length(optIdx)>1
            fprintf('Multiple points have same score function. Last point in tie list chosen.')
            optIdx = optIdx(end);
        end
        point = expDesigns(1,optIdx);
        pointType = expDesigns(2,optIdx);
        

       
        %Add point to current list of points for next round of calibration
        idx = find(fulldata.xdata==point);
        if (pointType == 1)
            currentPts = [currentPts; point fulldata.psadata(idx) pointType]
        else
            currentPts = [currentPts; point fulldata.anddata(idx) pointType]
        end
        
        %Remove skipped points from experimental design list
        beforePts = find(expDesigns(1,:)<point);
        expDesigns(:,beforePts) = [];
        
        %Remove chosen point
        chosenIdx = intersect(find(expDesigns(1,:)==point),find(expDesigns(2,:)==pointType));
        expDesigns(:,chosenIdx) = [];
        
        clear score relMI miVals unifMI
    end
end



point_list = currentPts; %Final list of points in order of selection


% Save results for later comparison
file = [path 'Results.mat'];
save(file, 'point_list','q2_vec','gamma1_vec','A0_vec','sigma2_vec','fixedVals','err_vec','pen_vec','penCoeffPSA_vec','penCoeffAND_vec','areaCI','fitChain')


end



%% Helper functions


% SSQ function for calibration of RT parameter
function SSrt = ssq_PCmodel_BothData(params,data,dayswitch,Ainit,Pinit)

global currentPts fixedVals

paramSet = [fixedVals(1:2) params(1) fixedVals(4:6) params(2) fixedVals(8) params(3) fixedVals(10:13) params(4) fixedVals(15:18)];

%modelout has androgen in col 1, PSA in col 2
[time,modelout] = evaluatePCmodel_BothData(paramSet,data,dayswitch,Ainit,Pinit);

%PSA set:
psaIdx = find(currentPts(:,3)==1);
psaoutput = interp1(time,modelout(:,1), currentPts(psaIdx,1));

%AND set:
andIdx = find(currentPts(:,3)==2);
andoutput = interp1(time,modelout(:,2), currentPts(andIdx,1));

SSrt = sum((psaoutput - currentPts(psaIdx,2)).^2)+sum((andoutput - currentPts(andIdx,2)).^2);

end

% Function for credible interval plotting
function v=predPCmodel_BothData(timef,params,data,dayswitch,Ainit,Pinit)

global fixedVals

paramSet = [fixedVals(1:2) params(1) fixedVals(4:6) params(2) fixedVals(8) params(3) fixedVals(10:13) params(4) fixedVals(15:18)];


[tsol, ysol] = evaluatePCmodel_BothData(paramSet,data,dayswitch,Ainit,Pinit);

v = interp1(tsol,ysol,timef);

end



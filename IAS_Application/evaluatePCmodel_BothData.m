function [tt,modelout] = evaluatePCmodel_BothData(params,data,dayswitch,Ainit,Pinit)

% Solve the ODE system
IC = [params(16) params(17)*params(16) params(18)*Ainit Ainit Pinit];
[tt,yy] = ode45(@(t,y)modelFunction(t,y,params,dayswitch),data.xdata,IC);

modelout = [yy(:,5) yy(:,4)]; %Return PSA in col 1, androgen in col 2
end
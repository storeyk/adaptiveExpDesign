function [tt,psaout] = evaluatePCmodel(params,data,dayswitch,Ainit,Pinit)

% Solve the ODE system
IC = [params(16) params(17)*params(16) params(18)*Ainit Ainit Pinit];
[tt,yy] = ode45(@(t,y)modelFunction(t,y,params,dayswitch),data.xdata,IC);

psaout = yy(:,5);
end
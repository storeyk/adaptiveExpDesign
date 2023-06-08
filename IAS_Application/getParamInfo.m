% Returns the ranges and mean values for parameters in the Meade model

function [lb,ub,meanVal,initGuess] = getParamInfo

mu = [.001 .09];
q1 = [.41 1.73];
q2 = [.01 .41];
d = [.001 .3];
c = [.00001 .0002]; 
K = [.8 1.7];
gamma1 = [.008 .8];
gamma2 = [.001 .1];
A0 = [2 22]; %pm 10 from max observed value
delta = [.03 .15];
mm = [.01 .9];
b = [.0001 .1];
sigma1 = [.001 1];
sigma2 = [.001 1];
epsilon = [.0001 .1];
x1init = [.009 .02];
x2initfac = [0 1];
Q0fac = [.5 1];

lb = [mu(1) q1(1) q2(1) d(1) c(1) K(1) gamma1(1) gamma2(1) A0(1) delta(1) mm(1) b(1) sigma1(1) sigma2(1) epsilon(1) x1init(1) x2initfac(1) Q0fac(1)];
ub = [mu(2) q1(2) q2(2) d(2) c(2) K(2) gamma1(2) gamma2(2) A0(2) delta(2) mm(2) b(2) sigma1(2) sigma2(2) epsilon(2) x1init(2) x2initfac(2) Q0fac(2)];
meanVal = lb+(ub-lb)/2;

%These give a good fit to PSA first 1.5 cycles 
initGuess = [.0154 .41 .01 .1 .00015 1 .3 .005 12 .1 .08 .02 .01 .01 .0001 .01 .0001 .8];

end
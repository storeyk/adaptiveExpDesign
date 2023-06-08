
function [tsol, ysol] = tumorModel_postRT(modelType,params)

%Specify initial condition
tsol = 0;
y = .02;
ysol = y;


treat = [0 15:19 22:26 29:33 36:40 43:47 50:54];

beta = params(end); %set alpha, estimate alpha/beta ratio
alpha = 1.5*beta;
d = 2; %dosage level
for nt = 1:(length(treat)-1)
    
    switch modelType 
        case 'EXP+DVR'
            [t,y] = ode23(@(t,y)tumorExponential(t,y,params),[treat(nt):1:treat(nt+1)],y(end));
        case {'LOG+DVR','LOG+PSI'}
            [t,y] = ode23(@(t,y)tumorLogistic(t,y,params),[treat(nt):1:treat(nt+1)],y(end));
        case 'LOG+CCR'
            [t,y] = ode23(@(t,y)tumorLogistic(t,y,params),[treat(nt):1:treat(nt+1)],y(end));
             params(2) = params(2)*(1-params(3)); 
    end
    ysol = [ysol; y(2:end)]; tsol = [tsol; t(2:end)];
    
    switch modelType 
        case {'EXP+DVR','LOG+DVR'}
            Radio = (1 - exp( (-alpha * d - beta * d^2 )))  * y(end);
            y(end) = y(end) - Radio;
        case {'LOG+PSI'} 
            K = params(2); 
            Radio = (1 - exp( (-alpha * d - beta * d^2 ))) * y(end) *(1 - y(end)/K);            
            y(end) = y(end) - Radio;
    end
    
end

end



function dval = tumorLogistic( t, val, params )

lambda = params(1);
K = params(2);

V = val(1);

dval(1) = lambda * V * ( 1 - V / K );

end



function dval = tumorExponential( t, val, params )

lambda = params(1);

V = val(1);

dval(1) = lambda * V;

end
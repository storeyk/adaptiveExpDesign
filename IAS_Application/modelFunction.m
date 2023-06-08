function ydot=modelFunction(t,y,params,dayswitch)
    %Model based on Meade, Weber, Phan et al. "High Accuracy Indicators
    % of Androgen Suppresion Therapy Failure for Prostate Cancer"

    
    mu = params(1);
    q1 = params(2);
    q2 = params(3);
    d = params(4);
    c = params(5);
    K = params(6);
    gamma1 = params(7);
    gamma2 = params(8);
    A0 = params(9);
    delta = params(10);
    mm = params(11);
    b = params(12);
    sigma1 = params(13);
    sigma2 = params(14);
    epsilon = params(15);
    

    x1=y(1); x2=y(2); Q=y(3); A=y(4); P=y(5);

    if mod(find(t<dayswitch,1),2)==0
        u=1;
    else
        u=0;
    end


    %Model Equations:
    dx1= max(mu*(1-q1/Q)*x1,0)-d*x1*(x1+x2)-c*K/(Q+K)*x1;
    dx2= max(mu*(1-q2/Q)*x2,0)-d*x2*(x1+x2)+c*K/(Q+K)*x1;
    dQ = mm*(A-Q)-(mu*(Q-q1)*x1+mu*(Q-q2)*x2)/(x1+x2);
    dA = gamma1*u*(1-A/A0)+gamma2-delta*A;
    dP = b*Q+max(sigma1*(1-q1/Q)*x1,0)+max(sigma2*(1-q2/Q)*x2,0)-epsilon*P;
    
    ydot = [dx1; dx2; dQ; dA; dP];
   
end


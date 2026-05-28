function dydt = HeterogeneousSingleLevelApprox_ODE(t,y,pars)
    
    beta = pars(1);
    theta = pars(2);
    eta = pars(3);
    gamma_i = pars(4);
    gamma_p = pars(5);
    N = pars(6);
    maxk = pars(7);
    pi_E_demon = pars(8);

    U = y(1:(maxk+1));
    E = y((maxk+2):(2*(maxk+1)));
    D = y((2*(maxk+1)+1):(3*(maxk+1)));
    P = y((3*(maxk+1)+1):(4*(maxk+1)));
    R = y((4*(maxk+1)+1):(5*(maxk+1)));

    pi_E = sum([0:maxk]'.*E)/pi_E_demon;

    dUdt = -beta*[0:maxk]'.*U*pi_E - theta*U.*P/N;
    dEdt = beta*[0:maxk]'.*U*pi_E + theta*U.*P/N - (eta+gamma_i)*E;
    dDdt = (eta+gamma_i)*E;
    dPdt = eta*E - gamma_p*P;
    dRdt = gamma_p*P;

    dydt = [dUdt;dEdt;dDdt;dPdt;dRdt];
end
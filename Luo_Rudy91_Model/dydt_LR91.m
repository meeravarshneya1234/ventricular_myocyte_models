function deriv = dydt_LR91(t,statevar,Id,p,c)

statevarcell = num2cell(statevar) ;
[V,m,h,j,d,f,n,Cai] = deal(statevarcell{:}) ;
 
  %%% compute ionic currents
  gNa = c.GNa_*m^3*h*j;
  INa = gNa*(V-p.ENa);

  Esi = 7.7 - 13.0287*log(Cai);
  gsi = c.Gsi_*d*f; 
  Isi = gsi*(V-Esi);

  if V > -100
    if ((V+77)<10^(-6))      % singularity
      ni = 0.7;
    else
      ni = 2.837*(exp(0.04*(V+77))-1)/((V+77)*exp(0.04*(V+35)));
    end
  else
    ni = 1;
  end
  gK = c.GK_*n*ni; 
  IK = gK*(V-p.EK);

  ax = 1.02/(1+exp(0.2385*(V-p.EK1-59.215)));
  bx = (0.49124*exp(0.08032*(V-p.EK1+5.476)) + ...
      exp(0.06175*(V-p.EK1-594.31)))/ ...
      (1+exp(-0.5143*(V-p.EK1+4.753)));
  x = ax/(ax+bx);
  gK1 = c.GK1_*x; 
  IK1 = gK1*(V-p.EK1);

  Kp = 1/(1+exp((7.488-V)/5.98));
  gKp = c.GKp_*Kp;
  IKp = gKp*(V-p.EKp);

  Ib = c.Gb_*(V-p.Eb);


  %% compute rate constants to update gates
  am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
  bm = 0.08*exp(-V/11); 

  if V >= -40
    ah = 0;
    aj = 0;
    bh = 1/(0.13*(1+exp((V+10.66)/(-11.1))));
    bj = 0.3*exp(-2.535*10^(-7)*V)/(1+exp(-0.1*(V+32)));
  else
    ah = 0.135*exp((80+V)/(-6.8));
    aj = (-1.2714*10^5*exp(0.2444*V)-3.474*10^(-5)*exp(-0.04391*V))...
        *(V+37.78)/(1+exp(0.311*(V+79.23)));
    bh = 3.56*exp(0.079*V)+3.1*10^5*exp(0.35*V);
    bj = 0.1212*exp(-0.01052*V)/(1+exp(-0.1378*(V+40.14)));
  end

  dmdt = am - (am+bm)*m ;
  dhdt = ah - (ah+bh)*h ;
  djdt = aj - (aj+bj)*j ;
  
  ad = 0.095*exp(-0.01*(V-5))/(1+exp(-0.072*(V-5)));
  bd = 0.07*exp(-0.017*(V+44))/(1+exp(0.05*(V+44))); 
  
  af = 0.012*exp(-0.008*(V+28))/(1+exp(0.15*(V+28)));
  bf = 0.0065*exp(-0.02*(V+30))/(1+exp(-0.2*(V+30))); 

  dddt = ad - (ad+bd)*d ;
  dfdt = af - (af+bf)*f ;
  
  an = 0.0005*exp(0.083*(V+50))/(1+exp(0.057*(V+50)));
  bn = 0.0013*exp(-0.06*(V+20))/(1+exp(-0.04*(V+20))); 
  dndt = an - (an+bn)*n ;
 
  dCaidt = -10^(-4)*Isi + 0.07*(10^(-4)-Cai);

  I = Id+INa+Isi+IK+IK1+IKp+Ib;
  dVdt = -I/p.Cm;

  deriv = [dVdt;dmdt;dhdt;djdt;dddt;dfdt;dndt;dCaidt] ;

return
  
  

% Luo-Rudy Phase 1 (1991) model
% generic ventricular cell 
% highly simplified model
%    
%    t   time variable
%    V   membrane potantial

colors = repmat('krgbmc',1,1000) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constants in the "p" structure ARE NOT varied in a population
% constants in the "c" structure ARE varied in a population

% % Physical constants
p.F = 96.5 ;                   % Faraday's constant, C/mmol
p.R = 8.314 ;                  % gas constant, J/K
p.T = 273+37 ;                 % absolute temperature, K 
p.RTF = p.R*p.T/p.F ;
p.Cm = 1 ;                     % membrane capacitance, uF/cm^2;

% % Standard ionic concentrations for ventricular cells
p.Ko = 5.4;                  % mmol/L
p.Ki = 145;                  % mmol/L
p.Nao = 140;                 % mmol/L
p.Nai = 10;                  % mmol/L
p.Cao = 1.8;                 % mmol/L

% computed quantities that do not change during simulation 
p.ENa = p.RTF*log(p.Nao/p.Nai);             % Nernst potential of Na, mV
c.GNa_ = 16;                          % mS/cm^2
c.GNa_ = 50 ;

c.Gsi_ = 0.09;

c.GK1_ = 0.6047*sqrt(p.Ko/5.4);
p.EK1 = p.RTF*log(p.Ko/p.Ki);

c.GK_ = 0.282*sqrt(p.Ko/5.4);         % mS/uF
p.PNa_K = 0.01833;                  % permability ratio of Na to K
p.EK = p.RTF*log((p.Ko+p.PNa_K*p.Nao)/(p.Ki+p.PNa_K*p.Nai));        % mV

c.GKp_ = 0.0183;
p.EKp = p.RTF*log(p.Ko/p.Ki);

c.Gb_ = 0.03921;
p.Eb = -59.87;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 2:  Define simulation, stimulus, and recording parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCL = 1000 ;  % Interval bewteen stimuli,[ms]
stim_delay = 100 ; % Time the first stimulus, [ms]
stim_dur = 2 ; % Stimulus duration
stim_amp = 20 ; % Stimulus amplitude 
nBeats = 100 ; % Number of beats to simulate 

stim_starts = stim_delay + PCL*(0:nBeats-1)  ;
stim_ends = stim_starts + stim_dur ;

% Create intervals for each beat 
simints = 3*nBeats ;
for i=1:nBeats
    intervals(3*i-2,:) = [PCL*(i-1),stim_starts(i)] ; %beginning 
    intervals(3*i-1,:) = [stim_starts(i),stim_ends(i)] ; %stimulus 
    intervals(3*i,:) = [stim_ends(i),PCL*i] ; % stimulus ends 
end
tend = nBeats*PCL ;              % end of simulation, ms
intervals(end,:) = [stim_ends(end),tend] ;

% Determine when to apply stim_amp or 0 amp  
Istim = zeros(simints,1) ;
stimindices = 3*(1:nBeats) - 1 ; % apply stimulus on second part of intervals
Istim(stimindices) = -stim_amp ; 

ssbefore = 1;% Run the model with no stimulus for 60 seconds. (only run if not running to steady state beforehand)
numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:  Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ic.V = -84 ;
ic.m = 0 ;
ic.h = 1 ;
ic.j = 1 ;
ic.d = 0 ;
ic.f = 1 ;
ic.n = 0 ;
ic.Cai = 0.12e-4 ; % mM

y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_LR91;

statevar_i = y0;
%%% let model rest for 60 seconds
if (ssbefore)
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    [~,posstatevars] = ode15s(odefcn,[0,60000],statevar_i,options,0,p,c) ;
    statevar_i = posstatevars(end,:) ;
end

%%% stimulate cell
if (nBeats > numbertokeep)
    for i=1:simints-3*numbertokeep
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
        statevar_i = posstatevars(end,:) ;
        t = post(end) ;
    end % for
    statevars = statevar_i ;
    for i=simints-3*numbertokeep+1:simints
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end % for
else
    t = 0 ;
    statevars = statevar_i ;
    for i=1:simints
        options = odeset('RelTol',1e-3,'AbsTol',1e-6);
        [post,posstatevars] = ode15s(odefcn,intervals(i,:),statevar_i,options,Istim(i),p,c) ;
        t = [t;post(2:end)] ;
        statevars = [statevars;posstatevars(2:end,:)] ;
        statevar_i = posstatevars(end,:) ;
    end % for
end % if

% uncomment if want to plot other state variables 
%[V,m,h,j,d,f,n,Cai] = deal(outputcell{:}) ; 

t = t - min(t) ;
V = statevars(:,V_ind);
state_variables = num2cell(statevars,1) ;
APD= find_APD(t,V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 5:  Plot AP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(t,V,'linewidth',2)
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('time (ms)')
ylabel('Voltage (mV)')
title(['Luo-Rudy 91 @ ' num2str(PCL) 'ms'])

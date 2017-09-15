% Hund et al model
% 2004
% dog ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constants
p.F = 96485 ;                  % Faraday's constant, C/mol
p.R = 8314 ;                   % gas constant, mJ/K
p.T = 273+25 ;                 % absolute temperature, K
p.RTF = p.R*p.T/p.F ;

% Cell geometry constants
p.Acap = 1.534e-4 ;            % cm^2
p.Vmyo = 25.84e-6 ;            % uL
p.Vss = 0.76e-6 ;              % uL
p.VJSR = 0.182e-6 ;            % uL
p.VNSR = 2.098e-6 ;            % uL

% Fixed ionic concentrations
% Initial conditions of others listed below
p.Ko = 5.4 ;                  % uM
p.Nao = 140 ;               % uM
p.Clo = 100 ;
p.Cao = 1.8 ;                 % uM

% % New global variables that will be defined in 'hund_ode.m
c.GNa_= 8.25 ;                          % mS/cm^2
c.GNaL_= 6.5e-3 ;
c.PCa_ = 2.43e-4 ;                       % cm/s
c.PCab = 1.995084e-7 ;                  % cm/s
p.gamma_Cao = 0.341 ;                   % dimensionless
p.gamma_Cai = 1 ;                       % dimensionless
p.dtaufCa_max = 10 ;                    % ms
c.GK1_ = 0.5 ;
c.GKr_ = 0.0138542 ;
c.GKs_ = 0.0248975 ;
p.pKNa = 0.01833 ;                  % relative permeability of IKs, Na to K
c.GKp_ = 2.76e-3 ;
c.Gto_ = 0.19 ;
p.Kmto2 = 0.1502 ;
c.GClb_ = 2.25e-4 ;
p.CTKCl_ = 7.0756e-6 ;
p.CTNaCl_ = 9.8443e-6 ;
c.PCl = 4e-7 ;
p.KmNa_NaK = 10 ;             % Half-saturation concentration of NaK pump (mM)
p.KmK_NaK = 1.5 ;             % Half-saturation concentration of NaK pump (mM)
c.INaK_ = 0.61875 ;           % Max. current through Na-K pump (uA/uF)
c.kNaCa = 4.5 ;
p.ksat = 0.27 ;
p.eta = 0.35 ;
p.KmNai = 12.3 ;
p.KmNao = 87.5 ;
p.KmCai = 0.0036 ;
p.KmCao = 1.3 ;
p.KmCa_allo = 1.25e-4 ;
c.Krel_ = 3000 ;
% This Km is same for CaMK effects on release, ICa gating, and uptake
p.KmCaMK = 0.15 ;
p.dtau_rel_max = 10 ;
c.IpCa_ = 0.0575 ; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
p.KmpCa = 0.5e-3 ; % Half-saturation concentration of sarcolemmal Ca pump (mM)
c.Vup = 0.004375 ;
p.Kmup = 0.00092 ;
p.CaNSR_max = 15.0 ;
p.tau_diff = 0.2 ;
p.tau_transfer = 120 ;
p.dKmPLBmax = 0.00017 ;
p.dJupmax = 0.75 ;

% % Buffers in cytosol
p.TRPNtot = 70e-3 ;
p.KmTRPN = 0.5e-3 ;
p.CMDNtot = 50e-3 ;
p.KmCMDN = 2.38e-3 ;

% % Buffers in JSR
p.CSQNtot = 10 ;
p.KmCSQN = 0.8 ;

% % Buffers in subspace
p.BSRtot = 0.047 ;
p.KmBSR = 0.00087 ;
p.BSLtot = 1.124 ;
p.KmBSL = 0.0087 ;
p.CaMK0 = 0.05 ;
p.KmCaM = 0.0015 ;

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

ssbefore = 1;% Run the model with no stimulus for 60 seconds. 
numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:  Set initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ic.V = -94.7 ;
ic.Cai = 0.0822e-3 ;
ic.Cass = 0.0822e-3 ;
ic.CaJSR = 1.25 ;
ic.CaNSR = 1.25 ;
ic.Nai = 9.71 ;
ic.Ki = 142.82 ;
ic.Cli = 19.53 ;
ic.m = 2.46e-4 ;
ic.h = 0.99869 ;
ic.j = 0.99887 ;
ic.hL = 0.99869 ;
ic.d = 1e-4 ;
ic.dp = 1 ;
ic.f = 0.983 ;
ic.f2 = 0.983 ;
ic.fCa = 0.942 ;
ic.fCa2 = 0.942 ;
ic.xKr = 0.229 ;
ic.xs1 = 1e-4 ;
ic.xs2 = 1e-4 ;
ic.ato = 3.742e-5 ;
ic.ito1 = 1 ;
ic.ito2 = 1 ;
ic.aa = 0 ;
ic.ro = 0 ;
ic.ri = 0 ;
ic.CAMK_trap = 0.001 ;
y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Hund;

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

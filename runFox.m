% Fox et al model
% AJP 282:H516-H530, 2002
% dog ventricular cell
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constants
p.F = 96.5 ;                   % Faraday constant, coulombs/mmol
p.R = 8.314 ;                  % gas constant, J/K
p.T = 273+25 ;                 % absolute temperature, K
p.RTF = p.R*p.T/p.F ;

% Cell geometry constants
p.Acap = 1.534e-4 ;            % cm^2
p.Cm = 1 ;                     % uF/cm^2
p.Vmyo = 25.84e-6 ;            % uL
p.VSR = 2e-6 ;                 % uL

% Fixed ionic concentrations
% Initial conditions of others listed below
p.Ko = 4000 ;                  % uM
p.Nao = 138000 ;               % uM
p.Cao = 2000 ;                 % uM

p.Ki = 148400 ;                % uM
p.Nai = 10000 ;                % uM

c.GNa_ = 12.8 ;                  % ms/uF
c.GK1_ = 2.8 ;                   % ms/uF
c.GKr_ = 0.0136 ;                % ms/uF
c.GKs_ = 0.0245 ;                % ms/uF
c.GKp_ = 0.002216 ;               % ms/uF
c.Gto_ = 0.23815 ;               % ms/uF
c.GNab_ = 0.0031 ;               % ms/uF
c.GCab_ = 3.842e-4 ;             % ms/uF

c.PCa_ = 2.26e-5 ;                % cm/ms
c.PCaK = 5.79e-7 ;               % cm/ms
c.Prel = 6 ;                     % ms-1
c.Pleak = 1e-6 ;                 % ms-1

c.INaK_ = 0.693 ;                % uA/uF
p.ICahalf = -0.265 ;             % uA/uF
c.IpCa_ = 0.05 ;                 % uA/uF
c.kNaCa = 1500 ;                 % uA/uF
c.Vup = 0.1 ;                    % uM/ms

% % values needed for calculation of NCX current
p.eta = 0.35 ;
p.ksat = 0.2 ;

% Several saturation constants, nomenclature gets tricky
% NCX
p.KmNa = 87500 ;                 % uM
p.KmCa = 1380 ;                  % uM
% Ca-inactivation of L-type current
p.KmfCa = 0.18 ;                 % uM
% IK1
p.KmK1 = 13000 ;                 % uM
% Na-K pump
p.KmNai = 10000 ;                % uM
p.KmKo = 1500 ;                  % uM
% sarcolemmal Ca pump
p.KmpCa = 0.05 ;                 % uM
% SERCA
p.Kmup = 0.32 ;                  % uM

% Buffering parameters
p.CMDNtot = 10 ;                 % uM
p.KmCMDN = 2 ;                   % uM
p.CSQNtot = 10000 ;              % uM
p.KmCSQN = 600 ;                 % uM

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
ic.Cai = 0.0472 ;
ic.CaSR = 320 ;
ic.f = 0.983 ;
ic.d = 1e-4 ;
ic.m = 2.46e-4 ;
ic.h = 0.99869 ;
ic.j = 0.99887 ;
ic.fCa = 0.942 ;
ic.xKr = 0.229 ;
ic.xKs = 1e-4 ;
ic.xto = 3.742e-5 ;
ic.yto = 1 ;
y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Fox;

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

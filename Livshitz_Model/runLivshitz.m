% Livshitz et al model
% 2009
% guinea pig ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical constants
p.F = 96485 ;                  % Faraday's constant, C/mol
p.R = 8314 ;                   % gas constant, mJ/K
p.T = 273+37 ;                 % absolute temperature, K
p.RTF = p.R*p.T/p.F ;
p.FRT = 1/p.RTF ;

% Cell Geometry
p.length_cell = 0.01;       % Length of the cell (cm)
p.radius = 0.0011;     % Radius of the cell (cm)
p.Vcell = 1000*pi*p.radius*p.radius*p.length_cell ;     %   3.801e-5 uL   % Cell volume (uL)
p.Ageo = 2*pi*p.radius*p.radius+2*pi*p.radius*p.length_cell ;  %   7.671e-5 cm^2    % Geometric membrane area (cm^2)
p.Acap = 2*p.Ageo ;             %   1.534e-4 cm^2    % Capacitive membrane area (cm^2)
p.Vmyo = p.Vcell*0.68 ;    % Myoplasm volume (uL)
p.Vmito = p.Vcell*0.24 ;  % Mitochondria volume (uL)
% data.vsr = Vcell*0.06;    % SR volume (uL)
p.VNSR = p.Vcell*0.0552 ;   % NSR volume (uL)
p.VJSR =  p.Vcell*0.0048 ;   % JSR volume (uL)
p.Vss = p.Vcell*0.02 ;

% Fixed ionic concentrations
% Initial conditions of others listed below
p.Ko = 4.5 ;                  % uM
p.Nao = 140 ;               % uM
p.Cao = 1.8 ;                 % uM

% % Na current
c.GNa_= 16 ;                          % mS/cm^2
c.GNab = 0.004 ;
% GNaL_= 6.5e-3 ;

% % Ca current
c.PCa_ = 5.4e-4 ;                        % cm/s
c.PCa_Na = 6.75e-7 ;                    % cm/s
c.PCa_K = 1.93e-7 ;                     % cm/s
c.PCab = 1.995084e-7 ;                  % cm/s
p.gamma_Cao = 0.341 ;                   % dimensionless
p.gamma_Cai = 1 ;                       % dimensionless
p.gamma_Nao = 0.75 ;                    % dimensionless
p.gamma_Nai = 0.75 ;                    % dimensionless
p.gamma_Ko = 0.75 ;                     % dimensionless
p.gamma_Ki = 0.75 ;                     % dimensionless
% hLca = 1 ;                            % dimensionless, Hill coefficient
% % Need to make sure this variable not used elsewhere
p.KmCa = 6e-4 ;                         % Half saturation constant, mM
% % T-type & background currents
c.GCaT = 0.05 ;
c.GCab = 0.003016 ;

c.GK1_ = 0.75;
c.GKr_ = 0.02614 ;
c.GKs_ = 0.433 ;
p.pKNa = 0.01833 ;                  % relative permeability of IKs, Na to K
c.GKp_ = 5.52e-3 ;

c.INaK_ = 2.25 ;             % Max. current through Na-K pump (uA/uF)
p.KmNa_NaK = 10 ;             % Half-saturation concentration of NaK pump (mM)
p.KmK_NaK = 1.5 ;             % Half-saturation concentration of NaK pump (mM)

c.kNCX = 0.00025 ;
p.ksat = 0.0001 ;
p.eta = 0.15 ;

c.alpha_rel = 0.125 ;
p.Krel_inf = 1 ;
p.hrel = 9 ;
p.beta_tau = 4.75 ;
p.Krel_tau = 0.0123 ;

c.IpCa_ = 1.15 ; % Max. Ca current through sarcolemmal Ca pump (uA/uF)
p.KmpCa = 5e-4 ; % Half-saturation concentration of sarcolemmal Ca pump (mM)

c.Vserca = 8.75e-3 ;              % mM/ms
p.Kmserca = 9.0e-4 ;              % mM
p.CaNSR_max = 15.0 ;
p.tau_transfer = 120 ;

% % Buffers in cytosol
p.TRPNtot = 70e-3 ;
p.KmTRPN = 0.5e-3 ;

% % % % troponin buffering might be more complicated now, need to check this
% trpnf = 40;   % forward  buffered in TRPN (mM)
% trpnb = 0.02;   % backward  TRPN (mM)
% cmdnf = 100;   % forward  buffered in TRPN (mM)
% cmdnb = 0.238;   % backward  TRPN (mM)

p.CMDNtot = 50e-3 ;
p.KmCMDN = 2.38e-3 ;

% % Buffers in JSR
p.CSQNtot = 10 ;
p.KmCSQN = 0.8 ;

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
ic.V = -84.7 ;
ic.Cai = 0.0822e-3 ;
ic.CaNSR = 1.25 ;
ic.CaJSR = 1.25 ;
ic.Nai = 9.71 ;
ic.Ki = 142.82 ;
ic.m = 2.46e-4 ;
ic.h = 0.99869 ;
ic.j = 0.99887 ;
ic.d = 1e-4 ;
ic.f = 0.983 ;
ic.b = 1e-4 ;
ic.g = 0.983 ;
ic.xKr = 0.229 ;
ic.xs1 = 1e-4 ;
ic.xs2 = 1e-4 ;
ic.Jrel = 1e-4 ;
y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Livshitz;

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

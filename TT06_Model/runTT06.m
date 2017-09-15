% TT06 et al model
% 2006
% human ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.celltype = 'endo';
% Physical constants
p.F = 96.4853415 ;              % Faraday's constant, coulombs/mmol
p.R = 8.314472 ;                % gas constant, J/(K mol)
p.T = 310.0;                    % absolute temperature, K
p.RTF = p.R*p.T/p.F ;

p.Cm = 2 ;                      % uF/cm2
p.Acap = 5.6297*3.280e-5 ;             % cm2
% % The factor of 5.6297 inserted as empirical correction factor
% % ten Tusscher-Panfilov paper is vague about volumes, capacitances, etc.
% % C code (from website) includes for instance, the variable
% % CAPACITANCE=0.185;
% % Units?  how is this used exactly?
% % To account for this, I computed change in [Ca] with unitary current
% % in both models.  Surface area needed to be increased by factor of 5.6297 to
% % make the two derivatives equal

p.Vmyo = 16.404 ;               % pL
p.VSR = 1.094 ;                 % pL
p.Vss = 0.05468 ;               % pL
% Paper says these in units of (um)^3, but this is obviously wrong
% 1 (um)^3 = 1 fL; myoplasmic volume more like 16 pL, or 16000 fL

p.Ko = 5400 ;                   % uM
p.Nao = 140000 ;                % uM
p.Cao = 2000 ;                  % uM

c.GNa_ = 14.838 ;               % nS/pF
c.PCa_ = 3.980e-5 ;             % cm uF-1 ms-1
c.GK1_ = 5.405 ;                % nS/pF
c.GKr_ = 0.153 ;                % nS/pF

p.GpK_ = 0.0146 ;               % nS/pF
c.GpCa_ = 0.1238 ;              % nS/pF
c.GNab_ = 2.9e-4 ;              % nS/pF
c.GCab_ = 5.92e-4 ;             % nS/pF

if strcmp(p.celltype,'endo')==1
    c.Gto_ = 0.073 ;                % nS/pF
    c.GKs_ = 0.392 ;                % nS/pF
elseif strcmp(p.celltype,'mid')==1
    c.Gto_ = 0.294 ;                % nS/pF
    c.GKs_ = 0.098 ;                % nS/pF
elseif strcmp(p.celltype,'epi')==1
    c.Gto_ = 0.294 ;                % nS/pF
    c.GKs_ = 0.392 ;                % nS/pF
else
    fprintf('Invalid cell type entered. Please re-enter cell type and try again.')
end

p.pKNa = 0.03 ;                  % relative permeability, Na to K

% Maximum rates of intracellular transport mechanisms
c.INaK_ = 2.724 ;                % pA/pF
c.kNaCa = 1000 ;                 % pA/pF
c.Iup_ = 6.375 ;                 % uM/ms

% % values needed for calculation of NCX current
% Several saturation constants, nomenclature gets tricky
% NCX
p.KmNa = 87500 ;                 % uM
p.KmCa = 1380 ;                  % uM
p.ksat = 0.1 ;                   % unitless
p.alpha_ncx = 2.5 ;              % unitless
p.eta = 0.35 ;                   % unitless, actually gamma in paper
% Na-K pump
p.KmNai = 40000 ;                % uM
p.KmKo = 1000 ;                  % uM
% Sarcolemmal Ca pump
p.KpCa = 0.5 ;                    % uM
% SERCA
p.Kmup = 0.25 ;                  % uM

% Vrel = 40.8 ;                 % ms^-1
c.Vrel = 0.102 ;                % ms^-1 as per KWTT source code
p.Vxfer = 0.0038 ;              % ms^-1
c.Vleak = 3.6e-4 ;              % ms^-1
% Note, units given in paper must be incorrect
% Paper says that these are in mM/ms
% But equations multiply rate by concentration gradient
% Thus no need to adjust rates for different concentration units

p.k1_ryr_prime = 0.15e-6 ;        % uM-2 ms-1
p.k2_ryr_prime = 0.045e-3 ;       % uM-1 ms-1
p.k3_ryr = 0.06 ;                 % ms-1
% k4_ryr = 1.5e-5 ;               % ms-1
p.k4_ryr = 0.005 ;                % ms-1 as per KWTT source code
p.maxsr = 2.5 ;                   % dimensionless
p.minsr = 1 ;                     % dimensionless
p.EC_ryr = 1500 ;                 % uM

% % % Buffering parameters
% % % Considered generic cytotolic and SR buffers
p.Bufc = 200 ;                  % uM
p.Kbufc = 1 ;                   % uM
p.Bufss = 400 ;                 % uM
p.Kbufss = 0.25 ;               % uM
p.BufSR = 10000 ;               % uM
p.KbufSR = 300 ;                % uM
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
ic.V = -83.5092 ;
ic.m = 0.0025 ;
ic.h = 0.6945 ;
ic.j = 0.6924 ;
ic.d = 4.2418e-005 ;
ic.f = 0.9697 ;
ic.f2 = 0.9784 ;
ic.fCass = 0.9999 ;
ic.r = 3.2195e-008 ;
ic.s = 1.0000 ;
ic.xs = 0.0038 ;
ic.xr1 = 3.1298e-004 ;
ic.xr2 = 0.4534 ;
ic.Rbar_ryr = 0.9816 ;
ic.Cai = 0.1061 ;
ic.Cass = 0.2381 ;
ic.CaSR = 3.6426e+003 ;
ic.Nai = 3.8067e+003 ;
ic.Ki = 1.2369e+005 ;

y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_TT06;

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

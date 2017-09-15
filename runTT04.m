% TT04 et al model
% 2004
% human ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.celltype = 'endo';
% Physical constants
p.F = 96.4867 ;                % Faraday constant, coulombs/mmol
p.R = 8.314 ;                  % gas constant, J/(K mol)
p.T = 273+37 ;                 % absolute temperature, K
p.RTF = p.R*p.T/p.F ;

p.Cm = 1 ;                      % uF/cm2
p.S = 0.2 ;                     % surface area to volume ratio

p.Vmyo = 16.404 ;               % pL
p.VSR = 1.094 ;                 % pL

p.Ko = 5400 ;                   % uM
p.Nao = 140000 ;                % uM
p.Cao = 2000 ;                  % uM

c.GNa_ = 14.838 ;               % nS/pF
c.PCa_ = 1.75e-4 ;              % cm3 uF-1 s-1
c.GK1_ = 5.405 ;                % nS/pF

c.GpK_ = 0.0146 ;               % nS/pF
c.GpCa_ = 0.025 ;               % nS/pF
c.GNab_ = 2.9e-4 ;              % nS/pF
c.GCab_ = 5.92e-4 ;             % nS/pF
c.GKr_ = 0.096 ;                % nS/pF
p.pKNa = 0.03 ;                  % relative permeability, Na to K

% Maximum rates of intracellular transport mechanisms
c.INaK_ = 1.362 ;                % pA/pF
c.kNaCa = 1000 ;                 % pA/pF
c.Iup_ = 0.425 ;                 % uM/ms
c.Vleak = 8e-5 ;                 % ms-1

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

c.arel = 16.464 ;                 % uM/ms
p.brel = 250 ;                   % uM
c.crel = 8.232 ;                  % uM/ms

% % % Buffering parameters
% % % Considered generic cytotolic and SR buffers
p.Bufc = 150 ;                % uM
p.Kbufc = 1 ;                   % uM
p.BufSR = 10000 ;              % uM
p.KbufSR = 300 ;                 % uM

if strcmp(p.celltype,'epi')==1
    % epicardial cell
    c.Gto_ = 0.294 ;                % nS/pF
    c.GKs_ = 0.245 ;                % nS/pF
elseif strcmp(p.celltype,'mid')==1
    % M cell
    c.Gto_ = 0.294 ;                % nS/pF
    c.GKs_ = 0.062 ;                 % nS/pF
elseif strcmp(p.celltype,'endo')==1
    % endocardial cell
    c.Gto_ = 0.073 ;                % nS/pF
    c.GKs_ = 0.245 ;                % nS/pF
else
    fprintf('Invalid cell type entered. Please re-enter cell type and try again.')
end
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
ic.V = -85 ;
ic.m = 0 ;
ic.h = 1 ;
ic.j = 1 ;
ic.d = 0 ;
ic.f = 1 ;
ic.fCa = 1 ;
ic.r = 0 ;
ic.s = 1 ;
ic.xs = 0 ;
ic.xr1 = 1 ;
ic.xr2 = 1 ;
ic.g = 0 ;
ic.Cai = 0.1 ;
ic.CaSR = 1000 ;
ic.Nai = 10000 ;
ic.Ki = 130000 ;

y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_TT04;

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

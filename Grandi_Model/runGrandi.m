% Grandi et al model
% 2010
% human ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.celltype = 'endo';
% Constants
p.R = 8314;       % [J/kmol*K]
p.Frdy = 96485;   % [C/mol]
p.Temp = 310;     % [K]
p.FoRT = p.Frdy/p.R/p.Temp;
p.Cmem = 1.3810e-10;   % [F] membrane capacitance
p.Qpow = (p.Temp-310)/10;

% Cell geometry
p.cellLength = 100;     % cell length [um]
p.cellRadius = 10.25;   % cell radius [um]
p.junctionLength = 160e-3;  % junc length [um]
p.junctionRadius = 15e-3;   % junc radius [um]
p.distSLcyto = 0.45;    % dist. SL to cytosol [um]
p.distJuncSL = 0.5;  % dist. junc to SL [um]
p.DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
p.DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
p.DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
p.DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec]
p.Vcell = pi*p.cellRadius^2*p.cellLength*1e-15;    % [L]
p.Vmyo = 0.65*p.Vcell; p.Vsr = 0.035*p.Vcell; p.Vsl = 0.02*p.Vcell; p.Vjunc = 0.0539*.01*p.Vcell;
p.SAjunc = 20150*pi*2*p.junctionLength*p.junctionRadius;  % [um^2]
p.SAsl = pi*2*p.cellRadius*p.cellLength;          % [um^2]
%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.3056e-11
% tau's from c-code, not used here
p.J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
p.J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
p.J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
p.J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
p.Fjunc = 0.11;   p.Fsl = 1-p.Fjunc;
p.Fjunc_CaL = 0.9; p.Fsl_CaL = 1-p.Fjunc_CaL;

% Fixed ion concentrations
p.Cli = 15;   % Intracellular Cl  [mM]
p.Clo = 150;  % Extracellular Cl  [mM]
p.Ko = 5.4;   % Extracellular K   [mM]
p.Nao = 140;  % Extracellular Na  [mM]
p.Cao = 1.8;  % Extracellular Ca  [mM]1.8
p.Mgi = 1;    % Intracellular Mg  [mM]

% Na transport parameters
p.KmNaip = 11;         % [mM]11
p.KmKo =1.5;         % [mM]1.5
p.Q10NaK = 1.63;
p.Q10KmNai = 1.39;

% Cl current parameters
p.KdClCa = 100e-3;    % [mM]

% I_Ca parameters
%%%% EAS: eventually modify this so that pNa/PCa_=constant
%%%%%     for now do not worry about it
p.pNa = 0.50*1.5e-8;       % [cm/sec]
p.pK = 0.50*2.7e-7;        % [cm/sec]
p.Q10CaL = 1.8;

% %% Ca transport parameters
p.KmCai = 3.59e-3;    % [mM]
p.KmCao = 1.3;        % [mM]
p.KmNai = 12.29;      % [mM]
p.KmNao = 87.5;       % [mM]
c.ksat = 0.32;        % [none]
p.nu = 0.27;          % [none]
p.Kdact = 0.150e-3;   % [mM]
p.Q10NCX = 1.57;      % [none]
p.KmPCa = 0.5e-3;     % [mM]
p.Q10SLCaP = 2.35;    % [none]

% SR flux parameters
p.Q10SRCaP = 2.6;          % [none]
p.Kmf = 0.246e-3;          % [mM] default
%Kmf = 0.175e-3;          % [mM]
p.Kmr = 1.7;               % [mM]L cytosol
p.hillSRCaP = 1.787;       % [mM]
p.koCa = 10;               % [mM^-2 1/ms]   %default 10   modified 20
p.kom = 0.06;              % [1/ms]
p.kiCa = 0.5;              % [1/mM/ms]
p.kim = 0.005;             % [1/ms]
p.ec50SR = 0.45;           % [mM]

% Buffering parameters
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
p.Bmax_Naj = 7.561;       % [mM] % Bmax_Naj = 3.7; (c-code difference?)  % Na buffering
p.Bmax_Nasl = 1.65;       % [mM]
p.koff_na = 1e-3;         % [1/ms]
p.kon_na = 0.1e-3;        % [1/mM/ms]
p.Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
p.koff_tncl = 19.6e-3;    % [1/ms]
p.kon_tncl = 32.7;        % [1/mM/ms]
p.Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity
p.koff_tnchca = 0.032e-3; % [1/ms]
p.kon_tnchca = 2.37;      % [1/mM/ms]
p.koff_tnchmg = 3.33e-3;  % [1/ms]
p.kon_tnchmg = 3e-3;      % [1/mM/ms]
p.Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
p.koff_cam = 238e-3;      % [1/ms]
p.kon_cam = 34;           % [1/mM/ms]
p.Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
p.koff_myoca = 0.46e-3;   % [1/ms]
p.kon_myoca = 13.8;       % [1/mM/ms]
p.koff_myomg = 0.057e-3;  % [1/ms]
p.kon_myomg = 0.0157;     % [1/mM/ms]
p.Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
p.koff_sr = 60e-3;        % [1/ms]
p.kon_sr = 100;           % [1/mM/ms]
p.Bmax_SLlowsl = 37.4e-3*p.Vmyo/p.Vsl;        % [mM]    % SL buffering
p.Bmax_SLlowj = 4.6e-3*p.Vmyo/p.Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
p.koff_sll = 1300e-3;     % [1/ms]
p.kon_sll = 100;          % [1/mM/ms]
p.Bmax_SLhighsl = 13.4e-3*p.Vmyo/p.Vsl;       % [mM]
p.Bmax_SLhighj = 1.65e-3*p.Vmyo/p.Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
p.koff_slh = 30e-3;       % [1/ms]
p.kon_slh = 100;          % [1/mM/ms]
p.Bmax_Csqn = 140e-3*p.Vmyo/p.Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
p.koff_csqn = 65;         % [1/ms]
p.kon_csqn = 100;         % [1/mM/ms]

c.GNa=23;
c.GNaB = 0.597e-3;    % [mS/uF] 0.897e-3
c.PCa_ = 0.50*5.4e-4;       % [cm/sec]
c.GCaB = 5.513e-4;    % [uA/uF] 3

c.GKr_ = 0.035 ;
c.GKs_ = 0.0035 ;
c.GK1_ = 0.35 ;
c.gkp = 2*0.001;
c.GClCa =0.5* 0.109625;   % [mS/uF]
c.GClB = 1*9e-3;        % [mS/uF]
c.IbarNaK = 1.0*1.8;%1.90719;     % [uA/uF]
c.IbarNCX = 1.0*4.5;      % [uA/uF]5.5 before - 9 in rabbit
c.ks = 25;                 % [1/ms]
c.Kleak_ = 5.348e-6 ;
c.Vmax_SRCaP = 1.0*5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
c.IbarSLCaP = 0.0673; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]

if strcmp(p.celltype,'epi')==1
    c.GtoSlow=1.0*0.13*0.12; %epi
    c.GtoFast=1.0*0.13*0.88; %epi0.88
elseif strcmp(p.celltype,'endo')==1
    c.GtoSlow=0.13*0.3*0.964; %endo
    c.GtoFast=0.13*0.3*0.036; %endo
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
ic.mo=1.405627e-3;
ic.ho= 9.867005e-1;
ic.jo=9.915620e-1;
ic.do=7.175662e-6;
ic.fo=1.000681;
ic.fcaBjo=2.421991e-2;
ic.fcaBslo=1.452605e-2;
ic.xtoso=4.051574e-3;
ic.ytoso=9.945511e-1;
ic.xtofo=4.051574e-3;
ic.ytofo= 9.945511e-1;
ic.xkro=8.641386e-3;
ic.xkso= 5.412034e-3;
ic.RyRro=8.884332e-1;
ic.RyRoo=8.156628e-7;
ic.RyRio=1.024274e-7;
ic.NaBjo=3.539892;
ic.NaBslo=7.720854e-1;
ic.TnCLo=8.773191e-3;
ic.TnCHco=1.078283e-1;
ic.TnCHmo=1.524002e-2;
ic.CaMo=2.911916e-4;
ic.Myoco=1.298754e-3;
ic.Myomo=1.381982e-1;
ic.SRBo=2.143165e-3;
ic.SLLjo=9.566355e-3;
ic.SLLslo=1.110363e-1;
ic.SLHjo=7.347888e-3;
ic.SLHslo=7.297378e-2;
ic.Csqnbo= 1.242988;
ic.Ca_sro=0.1e-1; %5.545201e-1;
ic.Najo=9.06;%8.80329;
ic.Naslo=9.06;%8.80733;
ic.Naio=9.06;%8.80853;
ic.Kio=120;
ic.Cajo=1.737475e-4;
ic.Caslo= 1.031812e-4;
ic.Caio=8.597401e-5;
ic.V=-8.09763e+1;
y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Grandi;

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

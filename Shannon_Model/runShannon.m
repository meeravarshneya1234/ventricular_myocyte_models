% Shannon et al model
% 2004
% rabbit ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Physical parameters
p.F = 96486.7;   % coulomb_per_mole (in model_parameters)
p.R2 = 8314.3;   % joule_per_kilomole_kelvin (R in model_parameters)
p.T = 310.0;   % kelvin (in model_parameters)

% % Cell geometry
p.Cm_per_area = 2.0e-6;   % farad_per_cm2 (in model_parameters)
p.cell_length = 100.0;   % micrometre (in model_parameters)
p.cell_radius = 10.25;   % micrometre (in model_parameters)
p.Fx_Ks_SL = 0.89;   % dimensionless (in IKs)
p.Fx_Ks_jct = 0.11;   % dimensionless (in IKs)
p.pKNa = 0.01833;   % dimensionless (in IKs)

% % Geometric and diffusional parameters
p.A_SL_cytosol = 1.3e-4;   % cm2 (in ion_diffusion)
p.A_jct_SL = 3.01e-6;   % cm2 (in ion_diffusion)
p.D_Ca_SL_cytosol = 1.22e-6;   % dm2_per_second (in ion_diffusion)
p.D_Ca_jct_SL = 1.64e-6;   % dm2_per_second (in ion_diffusion)
p.D_Na_SL_cytosol = 1.79e-5;   % dm2_per_second (in ion_diffusion)
p.D_Na_jct_SL = 1.09e-5;   % dm2_per_second (in ion_diffusion)
p.x_SL_cytosol = 0.45;   % micrometre (in ion_diffusion)
p.x_jct_SL = 0.5;   % micrometre (in ion_diffusion)

% % Ionic concentrations
p.Cao = 1.8;   % millimolar (in model_parameters)
p.Cli = 15.0;   % millimolar (in model_parameters)
p.Clo = 150.0;   % millimolar (in model_parameters)
p.Ki = 135.0;   % millimolar (in model_parameters)
p.Ko = 5.4;   % millimolar (in model_parameters)
p.Mgi = 1.0;   % millimolar (in model_parameters)
p.Nao = 140.0;   % millimolar (in model_parameters)

p.Fx_NaBk_SL = 0.89;   % dimensionless (in INab)
p.Fx_NaBk_jct = 0.11;   % dimensionless (in INab)
c.G_NaBk = 0.297e-3;   % milliS_per_microF (in INab)
p.Fx_Na_SL = 0.89;   % dimensionless (in INa)
p.Fx_Na_jct = 0.11;   % dimensionless (in INa)
c.G_INa = 16.0;   % milliS_per_microF (in INa)
c.G_tof = 0.02;   % milliS_per_microF (in Itof)
c.G_tos = 0.06;   % milliS_per_microF (in Itos)

p.Fx_Cl_SL = 0.89;   % dimensionless (in ICl_Ca)
p.Fx_Cl_jct = 0.11;   % dimensionless (in ICl_Ca)
c.G_Cl = 0.109625;   % milliS_per_microF (in ICl_Ca)
p.Kd_ClCa = 0.1;   % millimolar (in ICl_Ca)
c.G_ClBk = 0.009;   % milliS_per_microF (in IClb)

% % ICaL parameters
p.Fx_ICaL_SL = 0.1;   % dimensionless (in ICaL)
p.Fx_ICaL_jct = 0.9;   % dimensionless (in ICaL)
c.PCa_ = 5.4e-4;   % litre_per_farad_millisecond (in ICaL)
p.PK = 2.7e-7;   % litre_per_farad_millisecond (in ICaL)
p.PNa = 1.5e-8;   % litre_per_farad_millisecond (in ICaL)
p.Q10_CaL = 1.8;   % dimensionless (in ICaL)
p.gamma_Cai = 0.341;   % dimensionless (in ICaL)
p.gamma_Cao = 0.341;   % dimensionless (in ICaL)
p.gamma_Ki = 0.75;   % dimensionless (in ICaL)
p.gamma_Ko = 0.75;   % dimensionless (in ICaL)
p.gamma_Nai = 0.75;   % dimensionless (in ICaL)
p.gamma_Nao = 0.75;   % dimensionless (in ICaL)

c.GKr_ = 0.03 ;
c.GKs_ = 0.07 ;
c.GK1_ = 0.9 ;

% % SR Ca release parameters
c.KSRleak = 5.348e-6;   % per_millisecond (in Jleak_SR)
p.H2 = 1.787;   % dimensionless (H in Jpump_SR)
p.Kmf = 0.000246;   % millimolar (in Jpump_SR)
p.Kmr = 1.7;   % millimolar (in Jpump_SR)
p.Q10_SRCaP = 2.6;   % dimensionless (in Jpump_SR)
c.V_max_3 = 286.0e-6;   % millimolar_per_millisecond (V_max in Jpump_SR)
p.EC50_SR = 0.45;   % millimolar (in Jrel_SR)
p.HSR = 2.5;   % dimensionless (in Jrel_SR)
p.Max_SR = 15.0;   % dimensionless (in Jrel_SR)
p.Min_SR = 1.0;   % dimensionless (in Jrel_SR)
p.kiCa = 0.5;   % per_millimolar_per_millisecond (in Jrel_SR)
p.kim = 0.005;   % per_millisecond (in Jrel_SR)
p.koCa = 10.0;   % per_millimolar2_per_millisecond (in Jrel_SR)
p.kom = 0.06;   % per_millisecond (in Jrel_SR)
c.ks = 25.0;   % per_millisecond (in Jrel_SR)

% % Na-K pump paramters
p.Fx_NaK_SL = 0.89;   % dimensionless (in INaK)
p.Fx_NaK_jct = 0.11;   % dimensionless (in INaK)
p.H_NaK = 4.0;   % dimensionless (in INaK)
c.I_NaK_max = 1.91;   % microA_per_microF (in INaK)
p.Km_Ko = 1.5;   % millimolar (in INaK)
p.Km_Nai = 11.0;   % millimolar (in INaK)
p.Q10_Km_Nai = 1.49;   % dimensionless (in INaK)
p.Q10_NaK = 1.63;   % dimensionless (in INaK)

% % NCX parameters
p.Fx_NCX_SL = 0.89;   % dimensionless (in INaCa)
p.Fx_NCX_jct = 0.11;   % dimensionless (in INaCa)
p.HNa = 3.0;   % dimensionless (in INaCa)
p.K_mCai = 0.00359;   % millimolar (in INaCa)
p.K_mCao = 1.3;   % millimolar (in INaCa)
p.K_mNai = 12.29;   % millimolar (in INaCa)
p.K_mNao = 87.5;   % millimolar (in INaCa)
p.Kd_act = 0.000256;   % millimolar (in INaCa)
p.Q10_NCX = 1.57;   % dimensionless (in INaCa)
c.V_max_2 = 9.0;   % microA_per_microF (V_max in INaCa)
p.eta = 0.35;   % dimensionless (in INaCa)
c.ksat = 0.27;   % dimensionless (in INaCa)

% % SL Ca pump and background currents
p.Fx_CaBk_SL = 0.89;   % dimensionless (in ICab)
p.Fx_CaBk_jct = 0.11;   % dimensionless (in ICab)
c.G_CaBk = 0.0002513;   % milliS_per_microF (in ICab)
p.Fx_SLCaP_SL = 0.89;   % dimensionless (in ICap)
p.Fx_SLCaP_jct = 0.11;   % dimensionless (in ICap)
p.H1 = 1.6;   % dimensionless (H in ICap)
p.Km = 0.0005;   % millimolar (in ICap)
p.Q10_SLCaP = 2.35;   % dimensionless (in ICap)
c.V_maxAF = 0.0673;   % microA_per_microF (in ICap)

% % Buffering parameters
p.Bmax_Calsequestrin = 0.14;   % millimolar (in Ca_buffer)
p.Bmax_SLB_SL = 0.0374;   % millimolar (in Ca_buffer)
p.Bmax_SLB_jct = 0.0046;   % millimolar (in Ca_buffer)
p.Bmax_SLHigh_SL = 0.0134;   % millimolar (in Ca_buffer)
p.Bmax_SLHigh_jct = 0.00165;   % millimolar (in Ca_buffer)
p.koff_Calsequestrin = 65.0;   % per_millisecond (in Ca_buffer)
p.koff_SLB = 1.3;   % per_millisecond (in Ca_buffer)
p.koff_SLHigh = 30.0e-3;   % per_millisecond (in Ca_buffer)
p.kon_Calsequestrin = 100.0;   % per_millimolar_per_millisecond (in Ca_buffer)
p.kon_SL = 100.0;   % per_millimolar_per_millisecond (in Ca_buffer)

p.Bmax_SL = 1.65;   % millimolar (in Na_buffer)
p.Bmax_jct = 7.561;   % millimolar (in Na_buffer)
p.koff = 1.0e-3;   % per_millisecond (in Na_buffer)
p.kon = 0.0001;   % per_millimolar_per_millisecond (in Na_buffer)
p.Bmax_Calmodulin = 0.024;   % millimolar (in cytosolic_Ca_buffer)
p.Bmax_Myosin_Ca = 0.14;   % millimolar (in cytosolic_Ca_buffer)
p.Bmax_Myosin_Mg = 0.14;   % millimolar (in cytosolic_Ca_buffer)
p.Bmax_SRB = 0.0171;   % millimolar (in cytosolic_Ca_buffer)
p.Bmax_TroponinC = 0.07;   % millimolar (in cytosolic_Ca_buffer)
p.Bmax_TroponinC_Ca_Mg_Ca = 0.14;   % millimolar (in cytosolic_Ca_buffer)
p.Bmax_TroponinC_Ca_Mg_Mg = 0.14;   % millimolar (in cytosolic_Ca_buffer)
p.koff_Calmodulin = 238.0e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.koff_Myosin_Ca = 0.46e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.koff_Myosin_Mg = 0.057e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.koff_SRB = 60.0e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.koff_TroponinC = 19.6e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.koff_TroponinC_Ca_Mg_Ca = 0.032e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.koff_TroponinC_Ca_Mg_Mg = 3.33e-3;   % per_millisecond (in cytosolic_Ca_buffer)
p.kon_Calmodulin = 34.0;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
p.kon_Myosin_Ca = 13.8;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
p.kon_Myosin_Mg = 15.7e-3;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
p.kon_SRB = 100.0;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
p.kon_TroponinC = 32.7;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
p.kon_TroponinC_Ca_Mg_Ca = 2.37;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
p.kon_TroponinC_Ca_Mg_Mg = 3.0e-3;   % per_millimolar_per_millisecond (in cytosolic_Ca_buffer)
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
ic.Ca_Calsequestrin = 1.1865 ;
ic.Ca_SL = 0.0001064 ;
ic.Ca_SLB_SL = 0.0098686 ;
ic.Ca_SLB_jct = 0.0077808 ;
ic.Ca_SLHigh_SL = 0.11444 ;
ic.Ca_SLHigh_jct = 0.077504 ;
ic.Ca_SR = 0.54561 ;
ic.Ca_jct = 0.00017484 ;
ic.Cai = 8.735e-005 ;
ic.d = 6.9975e-006 ;
ic.fCaB_SL = 0.015353 ;
ic.fCaB_jct = 0.024609 ;
ic.f = 1.0007 ;
ic.Xr = 0.0084716 ;
ic.Xs = 0.006874 ;
ic.h = 0.98714 ;
ic.j = 0.99182 ;
ic.m = 0.0013707 ;
ic.X_tof = 0.0040112 ;
ic.Y_tof = 0.99463 ;
ic.R_tos = 0.38343 ;
ic.X_tos = 0.0040113 ;
ic.Y_tos = 0.29352 ;
ic.I = 9.272e-008 ;
ic.O = 7.1126e-007 ;
ic.R1 = 0.88467 ;
ic.Na_SL = 8.8741 ;
ic.Na_SL_buf = 0.77612 ;
ic.Na_jct = 8.8728 ;
ic.Na_jct_buf = 3.5571 ;
ic.Nai = 8.8745 ;
ic.V = -85.7197 ;
ic.Ca_Calmodulin = 0.00029596 ;
ic.Ca_Myosin = 0.0019847 ;
ic.Ca_SRB = 0.0021771 ;
ic.Ca_TroponinC = 0.0089637 ;
ic.Ca_TroponinC_Ca_Mg = 0.118 ;
ic.Mg_Myosin = 0.1375 ;
ic.Mg_TroponinC_Ca_Mg = 0.010338 ;

y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Shannon;

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

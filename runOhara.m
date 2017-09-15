% Ohara et al model
% 2011
% human ventricular cell 
%    
%    t   time variable
%    V   membrane potantial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1:  Define all constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.celltype = 'endo';
p.Nao=140.0;
p.Cao=1.8;
p.Ko=5.4;

%physical constants
p.R=8314.0;
p.T=310.0;
p.F=96485.0;
p.Cm=1.0;                     %uF

%cell geometry
p.L=0.01;
p.rad=0.0011;
p.vcell=1000*3.14*p.rad*p.rad*p.L;
p.Ageo=2*3.14*p.rad*p.rad+2*3.14*p.rad*p.L;
p.Acap=2*p.Ageo;
p.vmyo=0.68*p.vcell;
p.vnsr=0.0552*p.vcell;
p.vjsr=0.0048*p.vcell;
p.vss=0.02*p.vcell;
%
%CaMK constants
p.aCaMK=0.05;
p.bCaMK=0.00068;
p.CaMKo=0.05;

%reverse potential constant
p.PKNa=0.01833;

%INa constants
p.Ahf=0.99;
p.Ahs=1.0-p.Ahf;

%INaL constants
p.thL=200.0;

%Ito constants
p.delta_epi=1.0;
if strcmp(p.celltype,'epi')==1
    p.delta_epi=1.0-(0.95/(1.0+exp((V+70.0)/5.0)));
end

%ICaL, ICaNa, ICaK constants
p.Aff=0.6;
p.tjca=75.0;
p.Kmn=0.002;
p.k2n=1000.0;
p.zca=2.0;

%INaCa_i constants
p.kna1=15.0;
p.kna2=5.0;
p.kna3=88.12;
p.kasymm=12.5;
p.wna=6.0e4;
p.wca=6.0e4;
p.wnaca=5.0e3;
p.kcaon=1.5e6;
p.kcaoff=5.0e3;
p.qna=0.5224;
p.qca=0.1670;
p.zna=1.0;

%INaCa_ss constants
p.KmCaAct=150.0e-6;

%INaK constants
p.k1p=949.5;
p.k1m=182.4;
p.k2p=687.2;
p.k2m=39.4;
p.k3p=1899.0;
p.k3m=79300.0;
p.k4p=639.0;
p.k4m=40.0;
p.Knai0=9.073;
p.Knao0=27.78;
p.delta=-0.1550;

p.Kki=0.5;
p.Kko=0.3582;
p.MgADP=0.05;
p.MgATP=9.8;
p.Kmgatp=1.698e-7;
p.H=1.0e-7;
p.eP=4.2;
p.Khp=1.698e-7;
p.Knap=224.0;
p.Kxkur=292.0;

p.zk=1.0;

%calcium buffer constants
p.cmdnmax=0.05;
if strcmp(p.celltype,'epi')==1
    p.cmdnmax=p.cmdnmax*1.3;
end
p.kmcmdn=0.00238;
p.trpnmax=0.07;
p.kmtrpn=0.0005;
p.BSRmax=0.047;
p.KmBSR=0.00087;
p.BSLmax=1.124;
p.KmBSL=0.0087;
p.csqnmax=10.0;
p.kmcsqn=0.8;

%jsr constants
p.bt=4.75;
p.a_rel=0.5*p.bt;

% computed quantities that do not change during simulation
c.GNa=75;
c.GNaL=0.0075;
c.Gto=0.02;
c.GKr_=0.046;
c.GKs_=0.0034;
c.GK1=0.1908;
c.Gncx=0.0008;
c.GKb=0.003;
c.GpCa=0.0005;
c.PCa_=0.0001;
c.Pnak=30;

if  strcmp(p.celltype,'epi')==1
    c.GNaL=c.GNaL*0.6;
    c.Gto=c.Gto*4.0;
    c.GKr_=c.GKr_*1.3;
    c.GKs_=c.GKs_*1.4;
    c.GK1=c.GK1*1.2;
    c.Gncx=c.Gncx*1.1;
    c.GKb=c.GKb*0.6;
    c.PCa_=c.PCa_*1.2;
    c.Pnak=c.Pnak*0.9;
    
elseif  strcmp(p.celltype,'mid')==1
    c.Gto=c.Gto*4.0;
    c.GKr_=c.GKr_*0.8;
    c.GK1=c.GK1*1.3;
    c.Gncx=c.Gncx*1.4;
    c.PCa_=c.PCa_*2.5;
    c.Pnak=c.Pnak*0.7;
    
end
c.PNab=3.75e-10;
c.PCab=2.5e-8;

p.KmCaMK=0.15;
p.KmCaM=0.0015;

c.SERCA_total = 1 ;
c.RyR_total = 1 ;

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
ic.V=-87;
ic.Nai=7;
ic.Nass=ic.Nai;
ic.Ki=145;
ic.Kss=ic.Ki;
ic.Cai=1.0e-4;
ic.Cass=ic.Cai;
ic.Cansr=1.2;
ic.Cajsr=ic.Cansr;
ic.m=0;
ic.hf=1;
ic.hs=1;
ic.j=1;
ic.hsp=1;
ic.jp=1;
ic.mL=0;
ic.hL=1;
ic.hLp=1;
ic.a=0;
ic.iF=1;
ic.iS=1;
ic.ap=0;
ic.iFp=1;
ic.iSp=1;
ic.d=0;
ic.ff=1;
ic.fs=1;
ic.fcaf=1;
ic.fcas=1;
ic.jca=1;
ic.nca=0;
ic.ffp=1;
ic.fcafp=1;
ic.xrf=0;
ic.xrs=0;
ic.xs1=0;
ic.xs2=0;
ic.xk1=1;
ic.Jrelnp=0;
ic.Jrelp=0;
ic.CaMKt=0;
y0 = cell2mat(struct2cell(ic))';
V_ind=find(y0==ic.V); %determine where within state variables, the index of voltage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 4:  Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_Ohara;

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

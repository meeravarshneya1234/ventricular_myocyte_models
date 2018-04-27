%% Create a Virtual Population and Conduct PLS analysis 
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
%% Step 2:  Set up randomly varying parameter set  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = fieldnames(c); % names of parameters to vary 
n_parameters = length(F) ; % number of parameters to vary 
baselineparameters = c;
variations = 50 ; % number of model variants 
actionpotentials = cell(2*variations,1) ; % save voltage waveforms here 
sigma = 0.15 ;     % standard deviation of variation of parameters
allparameters = exp(sigma*randn(variations,n_parameters)) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 3:  Define simulation, stimulus, and recording parameters
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
%% Step 4:  Set initial conditions
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
%% Step 5: Run Simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
odefcn = @dydt_LR91;

for ii=1:variations
    scaling = allparameters(ii,:);        
    for iF = 1:length(F)
        aF = F{iF};
        c.(aF) = baselineparameters.(aF) * scaling(iF);
    end
    statevar_i = y0;
    
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
%     V = statevars(:,V_ind);
    outputcell = num2cell(statevars,1) ;
    [V,m,h,j,d,f,n,Cai] = deal(outputcell{:}) ;
    
    actionpotentials{2*ii-1} = t ;
    actionpotentials{2*ii} = V ;
    
    figure(gcf)
    plot(t,V,colors(ii))
    hold on
    grid on
    drawnow
    
    % save metrics for PLS analysis 
    APD(ii) = find_APD(t,V);
    Vrest(ii) = V(1) ;
    [peakV,~] = max(V) ;
    Vpeak(ii) = peakV ;    

    disp([int2str(ii),' of ',int2str(variations),' simulations completed.'])
    
%     datatable.APDs(ii,ind) = APD;
%     datatable.times{ii,ind} = t;
%     datatable.states{ii,ind} = statevars;
%     datatable.V{ii,ind} = V;
%     datatable.scalings(ii,:) = scaling;
    
end
figure(gcf)
xlabel('time')
ylabel('Voltage (mV)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 5: Run PLS  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = log(allparameters) ;
Y = [log(APD)',Vrest',Vpeak'] ;
logtransform = [1,0,0] ;

outputlabels = { ...
    'APD', ...
    'V_r_e_s_t', ...
    'V_p_e_a_k', ...
} ;

parameternames = {'G_N_a','G_s_i','G_K_1','G_K', ...
  'G_K_p','G_b', ...
} ;

% % Defining color table for intuitive display
igorct = zeros(256,3) ;
igorct(1,:) = 0 ;
igorct(256,:) = 0 ;
igorct(1:128,1) = 1 ;
igorct(1:128,2) = (1:128)/128 ;
igorct(1:128,3) = (1:128)/128 ;
igorct(129:256,1) = 1 - (1:128)/128 ;
igorct(129:256,2) = 1 - (1:128)/128 ;
igorct(129:256,3) = 1 ;

X_Z = zscore(X) ;
Y_Z = zscore(Y) ;
figure
handle = gcf ;
set(handle,'Position',[20,20,800,500])
set(handle,'PaperPosition',[1 1 8 5]) ;
subplot(1,2,1)
imagesc(X_Z)
colormap(igorct)
ylabel(['Model Variants n = ' num2str(variations)])
xticklabels(parameternames)
xtickangle(90)
subplot(1,2,2)
x = 1:1:length(outputlabels);
y = 1:variations; 
imagesc(x,y,Y_Z)
colormap(igorct)
xticks(1:1:length(outputlabels))
xticklabels(outputlabels)
xtickangle(90)
ylabel(['Model Variants n = ' num2str(variations)])
colorbar('Ticks',[-3 2],...
         'TickLabels',{'Decreased','Increased'})

nfact2keep = rank(X) ;
[T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y]=...
    PLS_nipals(X,Y,nfact2keep);

figure
handle = gcf ;
set(handle,'Position',[20,20,500,800])
set(handle,'PaperPosition',[1 1 5 8]) ;

for i=1:3

  subplot(3,1,i)
  hold on
  if logtransform(i)
    plot(exp(Y(:,i)),exp(Yhat(:,i)),'ro')
    plotmin = min(min([exp(Y(:,i)),exp(Yhat(:,i))])) ;
    plotmax = max(max([exp(Y(:,i)),exp(Yhat(:,i))])) ;
  else
    plot(Y(:,i),Yhat(:,i),'ro')
    plotmin = min([Y(:,i);Yhat(:,i)]) ;
    plotmax = max([Y(:,i);Yhat(:,i)]) ;
  end  
  xlabel(['True ',outputlabels{i}]) 
  ylabel(['Predicted ',outputlabels{i}])
  axis([plotmin plotmax plotmin plotmax])
  plot([plotmin,plotmax],[plotmin,plotmax],'k:')
  set(gca,'TickDir','Out') 

end

SSYT = sum((Y-ones(variations,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(variations,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

graycolor = [201/256*ones(256,1),201/256*ones(256,1),201/256*ones(256,1)] ;

figure
handle = gcf ;
set(handle,'Position',[20,20,500,800])
set(handle,'PaperPosition',[1 1 5 8]) ;

for i=1:3
  subplot(3,1,i)
  bar(Bpls(:,i))
  axis([0.5 n_parameters+0.5 -1 1]);
  xticklabels(parameternames)
  xtickangle(90)
  title(outputlabels{i})
  set(gca,'TickDir','out')
  colormap(graycolor)
end


%% Edit Settings 
settings.applyVoltageClamp = 0;
settings.vcRate = -100;
settings.runSignalingPathway = 1; 
settings.runElectrophysiol = 1; 
settings.APDRepLevel = 0.95; %used to calculate APD
settings.showProgress = 1; %shows progress depending on freq

settings.bcl = 1000; %ms
settings.freq = 3; %number of APs to run 
settings.LastBCL = settings.freq; 
settings.storeLast = NaN;  %stores all data for freq
settings.Istim = -80; 
settings.stimdur = 0.5;   
settings.ISO = 0;

settings.K_o = 5.4;
settings.Ca_o = 1.8; 
settings.Na_o = 140; 
settings.Cl_o = 100; 

settings.IKsNPParams = [7.3990e-003,3.1196e-002,8.0019e-001,5.6992e-003,4.1520e-002,1.3489e+000,3.8839e-001,-1.5019e-001,6.0693e-001,9.0654e-002,-1.1157e-001,2.8330e-002,3.1124e-003,-5.1660e-002,1.5522e+000,2.7304e-003,4.4198e-004,-1.2022e+000,4.0173e-004,2.0873e-004,1.9561e-001]; 
settings.IKsPParams =  [9.9415e-003,4.4809e-002,5.8172e-001,3.3201e-003,9.4217e-002,9.5364e-001,5.6356e-001,-1.7986e-001,5.8381e-001,6.5700e-002,-1.1899e-001,1.2406e-002,3.8525e-004,-6.4118e-002,7.7992e-001,4.6171e-003,2.3730e-004,-1.9742e+000,2.2652e-004,2.4689e-004,1.9561e-001];                                                                  
settings.LCCNPParams = [0.59,0.8,0.052,13,0.132,13.56,9.45,70,49.1,10.349,26.553,0.213,10.807,17.5,3,0.2474,13.825,6.3836,0.001,14.9186,1E-6,1.552E-4,1.1E-3,3.3696,1.2E-2]; 
settings.LCCPParams = [0.59,0.8,0.052,13,0.132,-4.7980,7.5699,70,49.1,10.349,38.494,0.213,10.807,29.979,3.1775,0.1,32.5,18.0,0.0001,6.0,1E-6,2.579E-04,0.002,10,0.01];
settings.INaNPParams = [2.15*8.25*1.1, 0.13, 10.66+16.7434, -11.1, 0.3, -2.535E-7, -0.1, 32, 0.135, 80+7, -6.8, 3.56, 0.079, 3.1E5, 0.35, -1.2714E5, 0.2444, -6.948E-5, -0.04391, 37.78, 0.311, 79.23, 0.1212, -0.01052, -0.1378, 40.14, -47.1-11.3729, -13.7299, -7]; 
settings.RyRP_Amp = 1.9925; 
settings.RyRP_Tau = 0.5357; 
settings.TauTr = 75; 
settings.observeICaLSS = 0; 

% Blocking Currents 
settings.ICaLB = 0.0; 
settings.IKsB = 0.0; 
settings.IKrB = 0.0; 
settings.INaKB =0.0; 
settings.INaCaB = 0.0;
settings.IKpB = 0.0; 
settings.IK1B = 0.0; 
settings.INabB = 0.0; 
settings.ITo1B = 0.0; 
settings.ITo2B = 0.0; 
settings.INaB = 0.0; 
settings.INaLB = 0.0; 
settings.IClB = 0.0; 
settings.IpCaB = 0.0; 
settings.CTKClB = 0.0; 
settings.CTNaClB = 0.0; 
settings.ICabB = 0.0; 
settings.IrelB = 0.0; 
settings.IupB = 0.0; 
settings.ItrB = 0.0; 
settings.IleakB = 0.0; 
settings.IdiffB = 0.0; 
settings.CAMKIIB = 0; 

%% Run Simulation 
[currents,State,Ti,APDs,settings]=mainHRdBA(settings); 

%% Plot all AP 
figure 
V = State(:,1);
plot(Ti,V,'linewidth',2)
set(gca,'FontSize',12,'FontWeight','bold')
xlabel('time (ms)')
ylabel('Voltage (mV)')

%% Plot all APs in freq individually 
% indc = find((mod(Ti,settings.bcl))==0);
% figure
% for i = 1:settings.freq
%     start = indc(2*i-1);
%     stop = indc(i*2);
%     
%     figure(gcf)
%     plot(Ti(start:stop),V(start:stop),'linewidth',2)
%     hold on 
% end 
% set(gca,'FontSize',12,'FontWeight','bold')
% xlabel('time (ms)')
% ylabel('Voltage (mV)')
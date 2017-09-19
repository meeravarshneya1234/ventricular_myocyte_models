% *************************************************************************
% *   Copyright (C) 2011 by Jordi Heijman, Paul G. A. Volders,            *
% *   Ronald L. Westra and Yoram Rudy                                     *
% *   Email: jordi.heijman@maastrichtuniversity.nl  / rudy@wustl.edu      *
% *   Web:   http://rudylab.wustl.edu                                     *
% *                                                                       *
% *   MATLAB Implementation of the model described in the article:        *
% *   "Local control of beta-adrenergic stimulation: Effects on           *
% *   ventricular myocyte electrophysiology and Ca2+ transient"           *
% *                                                                       *
% *   Journal of Molecular and Cellular Cardiology. 2011.                 *
% *                                                                       *
% *   Usage:                                                              *
% *   Settings can be stored in a structure and passed to the main        *
% *   function. These settings include (among other things):              *
% *     - bcl                 (pacing cycle length in ms)                 *
% *     - freq                (number of APs to calculate)                *
% *     - storeLast           (the number of beats that will be stored)   *
% *     - Istim               (amplitude of stimulus current)             *
% *     - stimdur             (duration of stimulus in ms)                *
% *     - X_o                 (extracellular concentration of ion X in mM)*
% *     - IXB                 (fraction of block of current X, e.g. ICaLB)*
% *     - ISO                 (concentration of isoproterenol in 10^-6 M) *
% *                                                                       *
% *   Example:                                                            *
% *     settings.bcl = 1000; settings.freq = 50; settings.storeLast = 5;  *
% *     settings.ISO = 1.0;                                               *
% *     [currents, State, t] = mainHRdBA(settings);                       *
% *     plot(t, State(:,1));                                              *
% *                                                                       *
% *   The program will return a structure with all individual currents, a *
% *   complete overview of the state vector over time and a list of the   *
% *   evaluated time points.                                              *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   funHRdBA.m, HRdBA_Currents.m                                        *
% *                                                                       *
% *   This program is free software; you can redistribute it and/or modify*
% *   it under the terms of the GNU General Public License as published by*
% *   the Free Software Foundation; either version 2 of the License, or   *
% *   (at your option) any later version.                                 *
% *                                                                       *
% *   This program is distributed in the hope that it will be useful,     *
% *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
% *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
% *   GNU General Public License for more details.                        *
% *                                                                       *
% *   You should have received a copy of the GNU General Public License   *
% *   along with this program; if not, write to the                       *
% *   Free Software Foundation, Inc.,                                     *
% *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.           *
% *************************************************************************

function [currents,State,Ti,APDs,settings]=mainHRdBA(settings)
    settings = setDefaultSettings(settings);
    settings.dataElectrophysiol = constants_Electrophysiol;
    settings.dataSignaling = constants_Signaling;
    
    if isfield(settings, 'applyMBC')
        if settings.applyMBC == 1
            settings.dataSignaling.J_CAV_ECAV = 1000 * settings.dataSignaling.J_CAV_ECAV;
            settings.dataSignaling.J_CAV_CYT = 1000 * settings.dataSignaling.J_CAV_CYT;
            settings.dataSignaling.J_ECAV_CYT = 1000 * settings.dataSignaling.J_ECAV_CYT; 
        end
    end
    if isfield(settings, 'applyPTX')
        if settings.applyPTX == 1, settings.dataSignaling.Gi_tot = 0.005; end
    end
    
    x0 = getInitialVector(settings); 
    State=[];Stm=[];Ti=[];
    currents = []; APDs = [];

    if isfield(settings, 'sparseJacobian') && ~isempty(settings.sparseJacobian)
        opts = odeset('RelTol', 1e-6,'JPattern','on');
    else
        opts= odeset('RelTol',1e-6, 'NonNegative', [2,3,4,5,6,13,80]);
    end
        
    if settings.applyVoltageClamp == 0
        % Current clamp mode: Loop through APs
        % ===========================================
        h4 = waitbar(0,' Matlab is calculating your data ...');
        for p=1:settings.freq                
            [t,X]=ode15s('funHRdBA',[0 settings.bcl],x0,opts,settings);

            if t(end) < 0.9 * settings.bcl
                disp(['Unable to complete AP, cancelling...']);
                break;
            end

            waitbar(p/settings.freq,h4)
            x0=[X(end,1:end)];

            % Record data
            if isnan(settings.storeLast) || p > settings.freq - settings.storeLast
                if isnan(settings.storeLast)
                    Ti=[Ti; t+(settings.bcl*(p-1))];
                else
                    Ti=[Ti; t+(settings.bcl*(p-1-(settings.freq - settings.storeLast)))];
                end
                State=[State; X];
                St=settings.Istim*ones(1,length(t));
                St((find(t>settings.stimdur)))=0;
                Stm=[Stm  St];
            end

            %curAPD = determineAPD(t, X(:,1), settings.APDRepLevel);
            curAPD = determineAPD(t, X(:,1));
            APDs = [APDs; curAPD];

            if settings.showProgress == 1
                disp(['Calculated AP ', num2str(p), ' / ', num2str(settings.freq), '... APD = ', num2str(curAPD)]);
            end
        end
        close(h4)        
    else
        % Voltage clamp mode
        % ===========================================
        if isfield(settings, 'voltageClampFile') && ~isempty(settings.voltageClampFile)
            settings.vcpFileData = dlmread(settings.voltageClampFile);            
        end
        
        if isfield(settings, 'vcpFileData') && ~isempty(settings.vcpFileData)
            settings.vcp = [0, settings.vcpFileData(end,1)];
        end
        
        prevt = 0;
        for p=1:size(settings.vcp,1)
            dur = settings.vcp(p,2) - prevt;        
            settings.Vhold = settings.vcp(p,1);
            
            if isfield(settings, 'applyJSRClamp') && p == settings.clampIndex
                x0(5) = settings.clampValue;
            end
        
            [t,X]=ode15s('funHRdBA',[0 dur],x0,opts,settings);

            if t(end) < 0.9 * dur
                disp(['Unable to complete VC, cancelling...']);
                break;
            end
            
            x0=[X(end,1:end)];
            Ti=[Ti; t(1:end-1)+(prevt)];
            State=[State; X(1:end-1, :)];
            prevt = prevt + dur;
        end
        Stm = zeros(length(Ti),1)';
    end
	currents = HRdBA_Currents(State,Ti,Stm,settings);     
end

% Sub-function to determine APD (between max dvdt and a given percentage of
% maximum and minimum membrane potential
% function APD = determineAPD(t, Vm, replevel)
%     dvdt = (Vm(2:end) - Vm(1:end-1)) ./ (t(2:end) - t(1:end-1));
%     [v, maxdvin] = max(dvdt);
%         
% 	baselinelvl = min(Vm);
% 	[peak, peakin] = max(Vm);
%     
% 	Voi = baselinelvl + (1 - replevel) * (peak - baselinelvl);
%     tin = find(Vm(peakin:end) < Voi);
%     if ~isempty(tin)
%         tin = tin(1);
%         APD = t(tin+peakin-1)-t(maxdvin);
%     else
%         APD = NaN;
%     end
% end
function APD = determineAPD(time, Voltage)
dVdt = diff(Voltage)./diff(time) ;
[~,dexmax] = max(dVdt) ;

%t of maximum dV/dt, consider this beginning of action potential
tinit = time(dexmax) ;

% Then determine peak V of action potential, for two reasons,
% 1) Because repolarization must, by definition, occur after this
% 2) To compute 50%, 90%, etc., must have this value
[~,peakdex] = max(Voltage) ;
tpeak = time(peakdex) ;
repoldex = find(time > tpeak & Voltage < -75,1) ; % APD based on repolarization back to -75 mV instead of -60 mV
if isempty(repoldex)
    APD = NaN;
else
    repoltime = time(repoldex(1)) ;
    APD = round((repoltime - tinit),1) ;
end

end

function Xinit = getInitialVector(settings)
      
        X0_Electro = [-8.7491E01, 1.3394E-02, 2.3413E-02, 2.3413E-02, 6.8659E+00, 1.1910E+00, 1.7546E-03, 6.8909E+00, 6.8909E+00, 1.4562E+02, 2.0273E+01, 2.0273E+01, 3.6675E-09, 6.8172E-03, 9.0163E-01, 9.9709E-01, 7.0530E-04, 3.6003E-01, 1.7687E-05, 9.9798E-01, 9.8747E-01, 1.2306E-08, 9.9604E-01, ...
              1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 9.1141E-01, 8.4395E-02, 2.9306E-03, 4.52285E-05, 2.6175E-07, 1.1424E-03, 7.9337E-05, 1.8366E-06, 1.4172E-08, 5.3696E-07, 2.4861E-08, 2.8776E-10, 1.1217E-10, 2.5967E-12, 8.7874E-15, 9.3722E-16, 1.6595E-17, 1.000, ...
              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 9.5624E-01, 4.2127E-02, 6.9597E-04, 5.1101E-06, 1.4070E-08, 9.0269E-04, 2.9826E-05, 3.2850E-07, 1.2060E-09, 3.1956E-07, 7.0390E-09, 3.8763E-11, 5.0277E-11, 5.5374E-13, 2.9664E-15, 1.1201E-15, 1.6613E-18, 1.2360E-3, 7.9472E-01, ...
              9.9123E-01, 6.8127E-04, 8.3805E-01, 9.9281E-01, 7.3074E-09, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]';
        
          X0_Signaling = zeros(57,1);
          X0_Signaling = [0.00685041638458665,0.0184627603007976,0.000731420577213056,0.00745773247314215,0.0191017987408719,0.00115141243826747,0.000607316088556676,0.000639038440072507,0.000419991861054322,0.347102959606005,9.62359241535767,0.474081735738211,0.0149041813757831,0.203016833596288,0.00944463350378086,2.49592854373432e-10,1.18055788874765e-09,7.07824478944671e-11,0.0904820284659604,0.00276490711096605,0.225475702283053,0.0326565916584703,0.192819110624505,0.205444874210056,0.174057375932567,0.817161796756964,0.567249910261073,0.249911886495890,0.0646928309115710,0.0664997605558791,0.489063888619456,0.362113356111496,0.126950532507959,0.0236821659448037,0.0128402905095188,0.00637363047239019,4.29171113639322e-05,0.00917039986149184,0.0282662056977524,0.000673713947839317,0.000765988420110534,0.592167467082831,0.673518785672382,0.239479458960528,0.126345311579566,0.00410693810508171,0.0589379755147718,0.0275455839709412,8.91799633266019e-10,0.00159632196638178,0.00209911481235842,0.000502792845976641,0.0110248953370551,1.13428924662652e-10,0.000364315164237569,0.000705656306851923,0.000341341142614041]';
        Xinit = [X0_Electro; X0_Signaling];
        
    if isfield(settings, 'X0')
        Xinit = settings.X0; 
    end
end
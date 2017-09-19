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
% *   This function will generate a selection of the BAR signaling        *
% *   validation incorporated into the model. It will create several dose *
% *   response curves for the major proteins / phosphorylation levels     *
% *                                                                       *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   mainHRdBA.m                                                   *
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

function showSignalingValidation
    close all;
    
    settings.bcl = 300000; % 5 minutes
    settings.freq = 6;
    settings.runElectrophysiol = 0;
    settings.storeLast = 6;
    
    ISO = logspace(-4,1,15)';
	[val, ind_500nM] = min(abs(ISO - 0.5));
    SS_Signaling = zeros(length(ISO),57);
    fICaLP = zeros(length(ISO),1); fIKsP = zeros(length(ISO),1);
    fPLBP = zeros(length(ISO),1); fTnIP = zeros(length(ISO),1);
    fINaP = zeros(length(ISO),1); fINaKP = zeros(length(ISO),1);
    fRyRP = zeros(length(ISO),1); fIKurP = zeros(length(ISO),1);
    cAMP = zeros(length(ISO),2); 
    State_500nM = []; t_500nM = []; currents_500nM = [];
    for i=1:length(ISO)
        settings.ISO = ISO(i);        
        [currents,State,Ti,APDs] = mainHRdBA(settings);
        
        SS_Signaling(i,:) = State(end,81:end);
        fICaLP(i) = currents.fICaLP(end); fIKsP(i) = currents.fIKsP(end);
        fPLBP(i) = currents.fPLBP(end); fTnIP(i) = currents.fTnIP(end);
        fINaP(i) = currents.fINaP(end); fINaKP(i) = currents.fINaKP(end);
        fRyRP(i) = currents.fRyRP(end); fIKurP(i) = currents.fIKurP(end);
        cAMP(i,:) = [0.02 * State(end,90) + 0.04 * State(end,91) + 0.678 * State(end,92), 0.02 * State(end,90) + 0.04 * State(end,91)];
        
        if i == ind_500nM
            currents_500nM = currents;
            State_500nM = State;
            t_500nM = Ti;
        end
    end

    ExpISOcAMP = [1.00E-06, 5.00E-04, 1.00E-03, 1.00E-02, 1.00E-01, 3.10E-01, 1.00E+00, 1.00E+01];
    ExpDatacAMPTotal = [0.330687831; 0.375661376; 0.380952381; 0.507936508; 0.80952381; 0.888888889; 0.970899471; 1];
    ExpDatacAMPParticulate = [0.132275132; 0.150793651; 0.174603175; 0.26984127; 0.35978836; 0.423280423; 0.44973545; 0.455026455];
    TOI_Combined = [0.1; 10; 20; 30; 40; 50; 60; 120; 300; 360; 600; 720; 1080];       
    Exp1_cAMP_Time = [0.314159292; 0.814159292; 1.0; 0.765486726; 0.612831858];
    Exp2_cAMP_Time = [0.281958296; 0.902085222; 1.0; 0.934723481; 0.912964642];
    
    figure;        
        subplot(1,2,1); semilogx(ISO, cAMP(:,1) ./ max(cAMP(:,1)), '-k', ISO, cAMP(:,2) ./ max(cAMP(:,1)), '-r', ExpISOcAMP, ExpDatacAMPTotal, 'ok', ExpISOcAMP, ExpDatacAMPParticulate, 'sr', 'LineWidth', 2, 'MarkerSize', 8); 
            title('Steady state cAMP'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Rel. cAMP Level');
            ylim([-0.1 1.1]);

        Rel_cAMP = 0.02 * State(:,90) + 0.04 * State(:,91) + 0.678 * State(:,92);
        subplot(1,2,2); plot(Ti / 1000, Rel_cAMP ./ max(Rel_cAMP), '-k', TOI_Combined([1; 8; 10; 12; 13]), Exp1_cAMP_Time, 'sr', TOI_Combined([1; 7; 8; 9; 11]), Exp2_cAMP_Time, 'ok', 'LineWidth', 2, 'MarkerSize', 8); 
            title('cAMP transient'); 
            xlabel('time (s)');
            ylabel('Rel. cAMP Level');
            xlim([-20 800]); ylim([-0.1 1.1]);
            
	ExpDataLCC = [0; 0.161752269; 0.647009077; 0.935350079];
    ExpISOLCC = [0,0.0100000000000000,0.100000000000000,1;];
    
    ExpDataIKS = [0.004085802; 0.116445352; 0.480081716; 0.877425945; 1];
    ExpISOIKS = [0,0.0100000000000000,0.0300000000000000,0.100000000000000,0.300000000000000;];
    
    ExpDataPLB = [0; 0.070986905; 0.429657582; 0.71920943; 0.926565915; 0.870523622; 1.04985896; 0.924697839; 0.909753227];
    ExpISOPLB = [0,0.00100000000000000,0.00500000000000000,0.0100000000000000,0.0250000000000000,0.0500000000000000,0.100000000000000,0.250000000000000,1;];
    
    ExpDataTnI = [0; 0.065076384; 0.11100185; 0.13075154; 0.288277627; 0.384878249; 0.580838504; 0.794086328; 0.907223599; 1; 0.944278262; 0.930500623; 0.951472672];
    ExpISOTnI = [0,0.00100000000000000,0.00250000000000000,0.00500000000000000,0.0100000000000000,0.0250000000000000,0.0500000000000000,0.0750000000000000,0.100000000000000,0.250000000000000,0.500000000000000,1,5;];
    
    ExpDataINa = [0; 0.09 / 0.3; 0.21 / 0.3; 0.29 / 0.3; 0.25 / 0.3];
    ExpISOINa = [0,0.00100000000000000,0.0100000000000000,0.100000000000000,1;];
    
    ExpDataINaK = [0; 0.128599385; 0.464076041; 0.948653434; 1.030658839];
    ExpISOINaK = [0,0.00100000000000000,0.0100000000000000,0.500000000000000,1;];
    
    ExpDataRyR = [0; 0.113855422; 0.514457831; 0.573493976; 0.877108434; 1.003614458];
    ExpISORyR = [0,0.0100000000000000,0.0500000000000000,0.100000000000000,0.300000000000000,1;];
    
    ExpDataIKur = [0; 0.095057034; 0.224334601; 0.623574144; 0.790874525; 1; 0.975285171];
    ExpISOIKur = [0,0.000100000000000000,0.00100000000000000,0.0100000000000000,0.0500000000000000,0.100000000000000,1];

    TOI = [0.01; 30; 60; 90; 120; 180; 240; 300];
    ExpDataICaLTime = [0; 0.1009; 0.6422; 0.8899; 0.9602; 0.9878];
    ExpDataPLBTime = [0; 0.473684211; 0.61754386; 0.835087719; 1; 0.694736842];
    ExpDataTnITime = [0; 0.253393665; 0.511312217; 0.814479638; 1; 0.959276018];

    figure;
        subplot(2,2,1); semilogx(ISO, fICaLP, '-k', ExpISOLCC, ExpDataLCC, 'ok', 'LineWidth', 2, 'MarkerSize', 8); 
            title('ICaL Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);
        
        subplot(2,2,2); semilogx(ISO, fRyRP, '-k', ExpISORyR, ExpDataRyR, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('RyR Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);
        
        subplot(2,2,3); semilogx(ISO, fINaP, '-k', ExpISOINa, ExpDataINa, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('INa Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);

        subplot(2,2,4); semilogx(ISO, fINaKP, '-k', ExpISOINaK, ExpDataINaK, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('INaK Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);

	figure;
        subplot(2,1,1); semilogx(ISO, fIKsP, '-k', ExpISOIKS, ExpDataIKS, 'ok', 'LineWidth', 2, 'MarkerSize', 8); 
            title('IKs Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);
        
        subplot(2,1,2); semilogx(ISO, fIKurP, '-k', ExpISOIKur, ExpDataIKur, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('IKur Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);
            
    figure;
        subplot(2,1,1); semilogx(ISO, fPLBP, '-k', ExpISOPLB, ExpDataPLB, 'ok', 'LineWidth', 2, 'MarkerSize', 8); 
            title('PLB Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);
        
        subplot(2,1,2); semilogx(ISO, fTnIP, '-k', ExpISOTnI, ExpDataTnI, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('TnI Phosphorylation'); 
            xlabel('[ISO] {\mu}mol/L');
            ylabel('Phosph. Level');
            ylim([-0.1 1.1]);
            
	figure;
        subplot(1,3,1); plot(t_500nM ./ 1000, currents_500nM.fICaLP, '-k', TOI([1,3,5,6,7,8]), ExpDataICaLTime, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('ICaL Phosphorylation Rate');
            xlabel('time (s)');
            ylabel('Phosph. Level');
            xlim([-10 300]); ylim([-0.1 1.1]);
            
        subplot(1,3,2); plot(t_500nM ./ 1000, currents_500nM.fPLBP, '-k', TOI([1,2,3,4,5,6]), ExpDataPLBTime, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('PLB Phosphorylation Rate');
            xlabel('time (s)');
            ylabel('Phosph. Level');
            xlim([-10 300]); ylim([-0.1 1.1]);
            
       subplot(1,3,3); plot(t_500nM ./ 1000, currents_500nM.fTnIP, '-k', TOI([1,2,3,4,5,6]), ExpDataTnITime, 'ok', 'LineWidth', 2, 'MarkerSize', 8);
            title('TnI Phosphorylation Rate');
            xlabel('time (s)');
            ylabel('Phosph. Level');
            xlim([-10 300]); ylim([-0.1 1.1]); 
end
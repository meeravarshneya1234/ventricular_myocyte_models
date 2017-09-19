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
% *   This function will generate a selection of the electrophysiological *
% *   validation incorporated into the model. It will create several I-V  *
% *   curves for the major currents under baseline conditions, ISO and/or *
% *   maximal CaMKII activation.                                          *
% *                                                                       *
% *   WARNING: This function may take a long time to compute!             *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   ComputeIVTraces.m                                                   *
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

function showElectroPhysiolValidation            
    
    %% ITo
    % ===================================
        outputfolder = '';
        settings = [];
        settings.runSignalingPathway = 0;
        settings.ISO = 0;
        
        settings.K_o = 4.0;
        settings.Ca_o = 1.8;
        
        settings.ICaLB = 1.0;
        settings.IKsB = 1.0;    
        settings.IKrB = 0.0;
        settings.INaKB = 1.0;
        settings.INaCaB = 1.0;
        settings.IKpB = 0.0;
        settings.IK1B = 0.0;
        settings.INabB = 1.0;
        settings.ITo1B = 0.0;
        settings.ITo2B = 1.0;
        settings.INaB = 1.0;
        settings.INaLB = 1.0;
        settings.IClB = 1.0;
        settings.IpCaB = 0.0;
        settings.CTKClB = 1.0;
        settings.CTNaClB = 1.0;
        settings.ICabB = 0.0;
        settings.IrelB = 1.0;
        settings.IupB = 1.0;
        settings.ItrB = 0.0;
        settings.IleakB = 1.0;
        settings.IdiffB = 0.0;
        
        % Voltage clamp protocols:
        % =================================================================           
        vcp_Recovery = [-80, 300; 40, 500; -80, 500; 40, 700; -80, 900];
        ind_Recovery = -3;
        val_Recovery = [10; 20; 30; 40; 50; 60; 75; 100; 150; 200; 300; 400; 500; 600; 800; 1000];
        rec_toi = [500, 700];

        settings.fIKsP = 0.0;
        settings.fPLBP = 0.0;
        settings.fRyRP = 0.0;
        settings.fICaLP = 0.0;
        settings.fINaP = 0.0;   
        settings.fINaKP = 0.0;
        settings.fTnIP = 0.0;
        settings.fIKurP = 0;        
        settings.fIToP_CaMKII = 0.0;
        
        [Traces_Recovery_Baseline, IV_ITo_Recovery_Baseline] = ComputeIVTraces(settings, vcp_Recovery, val_Recovery, ind_Recovery, 'ito', (000:1:3000)', rec_toi);        
        IV_ITo_Recovery_Baseline = IV_ITo_Recovery_Baseline' ./ IV_ITo_Recovery_Baseline(end);

        settings.fIKsP = 1.0;
        settings.fPLBP = 1.0;
        settings.fRyRP = 1.0;
        settings.fICaLP = 1.0;
        settings.fINaP = 1.0;   
        settings.fINaKP = 1.0;
        settings.fTnIP = 1.0;
        settings.fIKurP = 1.0;    
        settings.fIToP_CaMKII = 0.0;
        
        [Traces_Recovery_ISO, IV_ITo_Recovery_ISO] = ComputeIVTraces(settings, vcp_Recovery, val_Recovery, ind_Recovery, 'ito', (000:1:3000)', rec_toi);        
        IV_ITo_Recovery_ISO = IV_ITo_Recovery_ISO' ./ IV_ITo_Recovery_ISO(end);
        
        settings.fIKsP = 0.0;
        settings.fPLBP = 0.0;
        settings.fRyRP = 0.0;
        settings.fICaLP = 0.0;
        settings.fINaP = 0.0;   
        settings.fINaKP = 0.0;
        settings.fTnIP = 0.0;
        settings.fIKurP = 0;            
        settings.fIToP_CaMKII = 1.0;
        [Traces_Recovery_CaMKII, IV_ITo_Recovery_CaMKII] = ComputeIVTraces(settings, vcp_Recovery, val_Recovery, ind_Recovery, 'ito', (000:1:3000)', rec_toi);        
        IV_ITo_Recovery_CaMKII = IV_ITo_Recovery_CaMKII' ./ IV_ITo_Recovery_CaMKII(end);
            
        Exp_ITo_Recovery = [1.07E+01, 0.222423146;      1.96E+01, 0.215189873; ...
                            2.95E+01, 0.325497288;      4.11E+01, 0.428571429; ...
                            5.00E+01, 0.517179024;      5.98E+01, 0.591320072; ...
                            7.05E+01, 0.647377939;      8.04E+01, 0.665461121; ...
                            9.02E+01, 0.72875226;       1.00E+02, 0.732368897; ...
                            1.11E+02, 0.764918626;      1.21E+02, 0.763110307; ...
                            1.30E+02, 0.806509946;      1.41E+02, 0.808318264; ...
                            1.50E+02, 0.844484629;      1.56E+02, 0.839059675; ...
                            1.96E+02, 0.851717902;      2.36E+02, 0.87522604; ...
                            2.74E+02, 0.913200723;      3.16E+02, 0.893309222; ...
                            3.55E+02, 0.922242315;      3.95E+02, 0.938517179; ...
                            4.36E+02, 0.965641953;      4.74E+02, 0.956600362; ...
                            5.14E+02, 1.003616637;      5.55E+02, 1.001808318; ...
                            5.94E+02, 1.003616637;      6.35E+02, 1; ...
                            6.74E+02, 0.990958409];
 
        foptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-Inf; 0], 'StartPoint', [-1.0; 50]);
        ftype = fittype('a * exp(-x / b) + 1', 'options', foptions);
        [fitR,gofR,outR] = fit(val_Recovery, IV_ITo_Recovery_Baseline,ftype);        
        mdl_tau_ITo_Recovery_Baseline_SE = fitR.b;          
        [fitR,gofR,outR] = fit(val_Recovery, IV_ITo_Recovery_ISO,ftype);        
        mdl_tau_ITo_Recovery_ISO_SE = fitR.b; 
        [fitR,gofR,outR] = fit(val_Recovery, IV_ITo_Recovery_CaMKII,ftype);        
        mdl_tau_ITo_Recovery_CaMKII_SE = fitR.b;  
        
        foptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-Inf; 0; -Inf; 0], 'StartPoint', [-0.1; 10; -0.9; 500]);
        ftype = fittype('a * exp(-x / b) + c * exp(-x / d) + 1', 'options', foptions);
        [fitR,gofR,outR] = fit(val_Recovery, IV_ITo_Recovery_Baseline,ftype);        
        mdl_tau_ITo_Recovery_Baseline_DE = [fitR.b; fitR.d];        
        [fitR,gofR,outR] = fit(val_Recovery, IV_ITo_Recovery_ISO,ftype);        
        mdl_tau_ITo_Recovery_ISO_DE = [fitR.b; fitR.d];        
        [fitR,gofR,outR] = fit(val_Recovery, IV_ITo_Recovery_CaMKII,ftype);        
        mdl_tau_ITo_Recovery_CaMKII_DE = [fitR.b; fitR.d];
                    
        figure;
            subplot(1,3,1);
            plot(Exp_ITo_Recovery(:,1), Exp_ITo_Recovery(:,2), 'ok', val_Recovery, IV_ITo_Recovery_Baseline, '-k', val_Recovery, IV_ITo_Recovery_CaMKII, '-b', val_Recovery, IV_ITo_Recovery_ISO, '-r', 'LineWidth', 2, 'MarkerSize', 8);

            subplot(1,3,2);
            bar([mdl_tau_ITo_Recovery_Baseline_DE, mdl_tau_ITo_Recovery_ISO_DE, mdl_tau_ITo_Recovery_CaMKII_DE]);
            
            subplot(1,3,3);
            bar([mdl_tau_ITo_Recovery_CaMKII_DE(1) / mdl_tau_ITo_Recovery_Baseline_DE(1), 0.533; mdl_tau_ITo_Recovery_CaMKII_DE(2) / mdl_tau_ITo_Recovery_Baseline_DE(2), 0.344]);
            

        fid = fopen([outputfolder, 'Electrophysiology_Mdl_ITo_Recovery.txt'], 'w');
        fprintf(fid, 'variables=t \t IV_ITo_Recovery_Baseline \t IV_ITo_Recovery_ISO \t IV_ITo_Recovery_CaMKII\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f\n', [val_Recovery, IV_ITo_Recovery_Baseline, IV_ITo_Recovery_ISO, IV_ITo_Recovery_CaMKII]');
        fclose(fid);
        
        fid = fopen([outputfolder, 'Electrophysiology_Exp_ITo_Recovery.txt'], 'w');
        fprintf(fid, 'variables=t \t IV_ITo_Recovery_Baseline\n');
        fprintf(fid, '%8.6f \t %8.6f \n', [Exp_ITo_Recovery]');
        fclose(fid);
        
        fid = fopen([outputfolder, 'Electrophysiology_Bar_ITo_MdlTau.txt'], 'w');
        fprintf(fid, '%d \t %8.6f \t %8.6f \t %8.6f\n', [[1; 2; 3; 5; 6; 7], [diag([mdl_tau_ITo_Recovery_Baseline_DE(1), mdl_tau_ITo_Recovery_ISO_DE(1), mdl_tau_ITo_Recovery_CaMKII_DE(1)]); diag([mdl_tau_ITo_Recovery_Baseline_DE(2), mdl_tau_ITo_Recovery_ISO_DE(2), mdl_tau_ITo_Recovery_CaMKII_DE(2)])]]');
        fclose(fid);
        
        fid = fopen([outputfolder, 'Electrophysiology_Bar_ITo_ValTau.txt'], 'w');
        fprintf(fid, '%d \t %8.6f \t %8.6f \t %8.6f\n', [[1; 2; 4; 5], [mdl_tau_ITo_Recovery_CaMKII_DE(1) ./ mdl_tau_ITo_Recovery_Baseline_DE(1), 0, 0; 0, 0.533, 0.0822; mdl_tau_ITo_Recovery_CaMKII_DE(2) ./ mdl_tau_ITo_Recovery_Baseline_DE(2), 0, 0; 0, 0.344, 0.07642]]');
        fclose(fid);
        
    %% ICaL
    % ===================================
        outputfolder = '';
        settings = [];
        settings.runSignalingPathway = 0;
        settings.ISO = 0;
           
        settings.ICaLB = 0.0;
        settings.IKsB = 1.0;    
        settings.IKrB = 1.0;
        settings.INaKB =1.0;
        settings.INaCaB = 1.0;
        settings.IKpB = 1.0;
        settings.IK1B = 1.0;
        settings.INabB = 1.0;
        settings.ITo1B = 1.0;
        settings.ITo2B = 1.0;
        settings.INaB = 1.0;
        settings.INaLB = 1.0;
        settings.IClB = 1.0;
        settings.IpCaB = 0.0;
        settings.CTKClB = 1.0;
        settings.CTNaClB = 1.0;
        settings.ICabB = 0.0;
        settings.IrelB = 0.0;
        settings.IupB = 0.0;
        settings.ItrB = 0.0;
        settings.IleakB = 0.0;
        settings.IdiffB = 0.0;
    
        % Voltage clamp protocols:
        % =================================================================
        vcp_ACT_80 = [-80, 100; 70, 400; -80, 700];
        ind_ACT_80 = 2;
        val_ACT_80 = [-70; -60; -50; -40; -30; -20; -15; -10; -5; 0.0001; 5; 10; 15; 20; 30; 40; 50; 60; 70];
                    
        vcp_INACT  = [-80, 100; 70, 600; 10, 2100];
        ind_INACT  = 2;
        val_INACT  = [-60; -50; -40; -30; -20; -10; 0.00001; 10; 20];

        %vcp_Recovery = [-40, 200; 5, 600; -40, 600; 5, 1000; -40, 1200];
        vcp_Recovery = [-80, 300; 10, 400; -80, 400; 10, 500; -80, 700];
        ind_Recovery = -3;
        val_Recovery = [1:1:9 10:10:90 100:100:500]';
        %rec_toi = [600, 1000];
        rec_toi = [400, 500];
            
        settings.fIKsP = 0.0;
        settings.fPLBP = 0.0;
        settings.fRyRP = 0.0;
        settings.fICaLP = 0.0;
        settings.fINaP = 0.0;   
        settings.fINaKP = 0.0;
        settings.fTnIP = 0.0;
        settings.fIKurP = 0;        
        
        settings.Ca_o = 1.8;
        [Traces_NP, IV_ICaL_ACT_Baseline] = ComputeIVTraces(settings, vcp_ACT_80, val_ACT_80, ind_ACT_80, 'ical', (000:0.25:700)', [100, 300]);
        IV_ICaL_PEND_Baseline   = min((Traces_NP(1575:1599,:)));        

        settings.Ca_o = 2.5;
        [Traces, IV_ICaL_INACT_Baseline] = ComputeIVTraces(settings, vcp_INACT, val_INACT, ind_INACT, 'ical', (000:0.25:3100)', [600, 800]);
        IV_ICaL_INACT_Baseline = IV_ICaL_INACT_Baseline ./ IV_ICaL_INACT_Baseline(1);               
        
%         settings.IrelB = 1.0;
        settings.Ca_o = 2.0;
        [Traces_Recovery_Baseline, IV_ICaL_Recovery_Baseline] = ComputeIVTraces(settings, vcp_Recovery, val_Recovery, ind_Recovery, 'ical', (000:0.5:2000)', rec_toi);        
        IV_ICaL_Recovery_Baseline = IV_ICaL_Recovery_Baseline' ./ IV_ICaL_Recovery_Baseline(end);
        
        Trace_Show_Baseline = ComputeIVTraces(settings, [-90,300;0,500], 0, 2, 'ical', (000:0.5:500)', [300, 500]);
        [v,in] = min(Trace_Show_Baseline(:,2));
        t_adj = Trace_Show_Baseline(in:end,1) - Trace_Show_Baseline(in);
        ICaL_adj = (Trace_Show_Baseline(in:end,2) - Trace_Show_Baseline(end,2)) ./ (Trace_Show_Baseline(end,2) - Trace_Show_Baseline(in,2));
        foptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0], 'StartPoint', [10]);
        ftype = fittype('-exp(-x / b)', 'options', foptions);
        [fitR,gofR,outR] = fit(t_adj, ICaL_adj,ftype);        
        mdl_tau_ICaL_Baseline = fitR.b;
        [v,in] = min(abs(Trace_Show_Baseline - 290));
        Trace_Show_Baseline = [Trace_Show_Baseline(in:end,1) - 300, Trace_Show_Baseline(in:end,2)];
        
        settings.fICaLP_CaMKII = 1.0;
        settings.Ca_o = 1.8;
        [Traces_NP, IV_ICaL_ACT_CaMKII] = ComputeIVTraces(settings, vcp_ACT_80, val_ACT_80, ind_ACT_80, 'ical', (000:0.25:700)', [100, 300]);
        IV_ICaL_PEND_CaMKII   = min((Traces_NP(1575:1599,:)));        

        settings.Ca_o = 2.5;
        [Traces, IV_ICaL_INACT_CaMKII] = ComputeIVTraces(settings, vcp_INACT, val_INACT, ind_INACT, 'ical', (000:0.25:3100)', [600, 800]);
        IV_ICaL_INACT_CaMKII = IV_ICaL_INACT_CaMKII ./ IV_ICaL_INACT_CaMKII(1); 
        
        settings.Ca_o = 2.0;
        [Traces_Recovery_CaMKII, IV_ICaL_Recovery_CaMKII] = ComputeIVTraces(settings, vcp_Recovery, val_Recovery, ind_Recovery, 'ical', (000:0.5:2000)', rec_toi);
        IV_ICaL_Recovery_CaMKII = IV_ICaL_Recovery_CaMKII' ./ IV_ICaL_Recovery_CaMKII(end);
        
        Trace_Show_CaMKII = ComputeIVTraces(settings, [-90,300;0,500], 0, 2, 'ical', (000:0.5:500)', [300, 500]);
        settings = rmfield(settings, 'fICaLP_CaMKII');
        [v,in] = min(Trace_Show_CaMKII(:,2));
        t_adj = Trace_Show_CaMKII(in:end,1) - Trace_Show_CaMKII(in);
        ICaL_adj = (Trace_Show_CaMKII(in:end,2) - Trace_Show_CaMKII(end,2)) ./ (Trace_Show_CaMKII(end,2) - Trace_Show_CaMKII(in,2));
        [fitR,gofR,outR] = fit(t_adj, ICaL_adj,ftype);        
        mdl_tau_ICaL_CaMKII = fitR.b;
        [v,in] = min(abs(Trace_Show_CaMKII - 290));
        Trace_Show_CaMKII = [Trace_Show_CaMKII(in:end,1) - 300, Trace_Show_CaMKII(in:end,2)];
        
        settings.IrelB = 1.0;
        settings.fIKsP = 1.0;
        settings.fPLBP = 1.0;
        settings.fRyRP = 1.0;
        settings.fICaLP = 1.0;
        settings.fINaP = 1.0;
        settings.fINaKP = 1.0;
        settings.fTnIP = 1.0;
        settings.fIKurP = 1.0;
        
        settings.Ca_o = 1.8;
        [Traces_P, IV_ICaL_ACT_ISO] = ComputeIVTraces(settings, vcp_ACT_80, val_ACT_80, ind_ACT_80, 'ical', (000:0.25:700)', [100, 300]);
        IV_ICaL_PEND_ISO   = min((Traces_P(1575:1599,:)));
        
        settings.Ca_o = 2.5;
        [Traces, IV_ICaL_INACT_ISO] = ComputeIVTraces(settings, vcp_INACT, val_INACT, ind_INACT, 'ical', (000:0.25:3100)', [600, 800]);
        IV_ICaL_INACT_ISO = IV_ICaL_INACT_ISO ./ IV_ICaL_INACT_ISO(1);
        
        settings.IrelB = 0.0;
        settings.Ca_o = 2.0;
        [Traces_Recovery_ISO, IV_ICaL_Recovery_ISO] = ComputeIVTraces(settings, vcp_Recovery, val_Recovery, ind_Recovery, 'ical', (000:0.5:2000)', rec_toi);
        IV_ICaL_Recovery_ISO = IV_ICaL_Recovery_ISO' ./ IV_ICaL_Recovery_ISO(end);

        Trace_Show_ISO = ComputeIVTraces(settings, [-90,300;0,500], 0, 2, 'ical', (000:0.5:500)', [300, 500]);
        [v,in] = min(Trace_Show_ISO(:,2));
        t_adj = Trace_Show_ISO(in:end,1) - Trace_Show_ISO(in);
        ICaL_adj = (Trace_Show_ISO(in:end,2) - Trace_Show_ISO(end,2)) ./ (Trace_Show_ISO(end,2) - Trace_Show_ISO(in,2));
        [fitR,gofR,outR] = fit(t_adj, ICaL_adj,ftype);        
        mdl_tau_ICaL_ISO = fitR.b;
        [v,in] = min(abs(Trace_Show_ISO - 290));
        Trace_Show_ISO = [Trace_Show_ISO(in:end,1) - 300, Trace_Show_ISO(in:end,2)];
        
        res.V_ICaL_ACT = val_ACT_80;
        res.V_ICaL_INACT = val_INACT;
        res.t_ICaL_Recovery = val_Recovery;
        
        res.IV_ICaL_ACT_Baseline = IV_ICaL_ACT_Baseline';
        res.IV_ICaL_PEND_Baseline = IV_ICaL_PEND_Baseline';
        res.IV_ICaL_INACT_Baseline = IV_ICaL_INACT_Baseline';
        res.IV_ICaL_Recovery_Baseline = IV_ICaL_Recovery_Baseline;
        
        res.IV_ICaL_ACT_CaMKII = IV_ICaL_ACT_CaMKII';
        res.IV_ICaL_PEND_CaMKII = IV_ICaL_PEND_CaMKII';
        res.IV_ICaL_INACT_CaMKII = IV_ICaL_INACT_CaMKII';        
        res.IV_ICaL_Recovery_CaMKII = IV_ICaL_Recovery_CaMKII;
        
        res.IV_ICaL_ACT_ISO = IV_ICaL_ACT_ISO';
        res.IV_ICaL_PEND_ISO = IV_ICaL_PEND_ISO';
        res.IV_ICaL_INACT_ISO = IV_ICaL_INACT_ISO';
        res.IV_ICaL_Recovery_ISO = IV_ICaL_Recovery_ISO;
        
        res.Exp_ICaL_IV_Baseline = [(-30:10:70)', [-5.705000162E-001; -8.392000198E-001; -1.852100015E+000; -3.735899925E+000; -5.114099979E+000; ...
            -5.172100067E+000; -4.064199924E+000; -2.844099998E+000; -1.820600033E+000; -1.077999949E+000; -4.758000076E-001]];
        
        res.Exp_ICaL_Inact_Baseline = [-40.6384, 1.0134; -35.8896, 1.0191; -30.8418, 1.0076; -26.3087, 0.9561; -21.1328, 0.8302; -15.9804, 0.4752; ...
            -10.3838, 0.2233; -5.4748, 0.0859; -0.1153, 0.0458; 4.373, 0.0344; 9.6877, 0.0344; 14.2081, 5.73E-03; 19.4971, 0.0172];       
        
        res.Exp_ICaL_Recovery_Baseline = [4.71E+00, 0.127208481, 0.003533569; ...
                                 2.47E+01, 0.4204947, 0.03180212; ...
                                 4.59E+01, 0.51590106, 0.045936396; ...
                                 6.47E+01, 0.586572438, 0.053003534; ...
                                 8.59E+01, 0.614840989, 0.053003534; ...
                                 1.06E+02, 0.657243816, 0.053003534; ...
                                 1.26E+02, 0.696113074, 0.056537102; ...
                                 1.45E+02, 0.706713781, 0.045936396; ...
                                 1.65E+02, 0.731448763, 0.035335689; ...
                                 1.85E+02, 0.731448763, 0.03180212; ...
                                 2.06E+02, 0.763250883, 0.042402827; ...
                                 2.26E+02, 0.791519435, 0.038869258; ...
                                 2.46E+02, 0.784452297, 0.028268551; ...
                                 2.66E+02, 0.812720848, 0.038869258; ...
                                 2.86E+02, 0.851590106, 0.042402827; ...
                                 3.06E+02, 0.840989399, 0.024734982; ...
                                 3.25E+02, 0.855123675, 0.028268551; ...
                                 3.45E+02, 0.865724382, 0.024734982; ...
                                 3.66E+02, 0.879858657, 0.042402827; ...
                                 3.86E+02, 0.890459364, 0.035335689];

        figure;
            subplot(2,6,[1,2]); plot(res.Exp_ICaL_IV_Baseline(:,1), res.Exp_ICaL_IV_Baseline(:,2), 'ok', res.V_ICaL_ACT, res.IV_ICaL_ACT_Baseline, '-k', res.V_ICaL_ACT, res.IV_ICaL_ACT_ISO, '-r', res.V_ICaL_ACT, res.IV_ICaL_ACT_CaMKII, '-b', 'LineWidth', 2, 'MarkerSize', 8);
                title('ICaL I-V'); ylabel('ICaL (pA/pF)'); xlabel('Vm (mV)');
                
            subplot(2,6,[3,4]); plot(res.Exp_ICaL_Inact_Baseline(:,1), res.Exp_ICaL_Inact_Baseline(:,2), 'ok', res.V_ICaL_INACT, res.IV_ICaL_INACT_Baseline, '-k', res.V_ICaL_INACT, res.IV_ICaL_INACT_ISO, '-r', res.V_ICaL_INACT, res.IV_ICaL_INACT_CaMKII, '-b', 'LineWidth', 2, 'MarkerSize', 8);
                title('ICaL Inactivation'); ylabel('Rel. ICaL'); xlabel('Vm (mV)');
                
        foptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0; -Inf; 0], 'StartPoint', [0.2; -20; 4]);
        ftype = fittype('(1-a)/(1 + exp((x - b)/c)) + a', 'options', foptions);
        [fitR,gofR,outR] = fit(res.V_ICaL_INACT,res.IV_ICaL_INACT_Baseline,ftype);        
        mdl_vhalf_base = fitR.b

        foptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [0; -Inf; 0], 'StartPoint', [0.2; -20; 4]);
        ftype = fittype('(1-a)/(1 + exp((x - b)/c)) + a', 'options', foptions);
        [fitR,gofR,outR] = fit(res.V_ICaL_INACT,res.IV_ICaL_INACT_ISO,ftype);        
        mdl_vhalf_iso = fitR.b;

        [mv_iso, ind_iso] = max(abs(res.IV_ICaL_ACT_ISO));
        [mv_base, ind_base] = max(abs(res.IV_ICaL_ACT_Baseline))
        [mv_camk, ind_camk] = max(abs(res.IV_ICaL_ACT_CaMKII))

        subplot(2,6,5); bar([mv_iso ./ mv_base, 4.1; mv_camk ./ mv_base, 1.3]);
            title('ICaL modulation'); ylabel('Fold incr. from baseline');
            set(gca,'XTick',1:2); set(gca,'XTickLabel',{'ISO','CaMKII'});
            
        subplot(2,6,6); bar([abs(res.V_ICaL_ACT(ind_iso) - res.V_ICaL_ACT(ind_base)), 15; abs(mdl_vhalf_iso - mdl_vhalf_base), 8.0]);
            title(' '); ylabel('V-shift w. ISO');
            set(gca,'XTick',1:2); set(gca,'XTickLabel',{'Peak I-V','V_{1/2} Inact'});
            
%         subplot(2,6,[7,8]); plot(res.Exp_ICaL_Recovery_Baseline(:,1), res.Exp_ICaL_Recovery_Baseline(:,2), 'ok', res.t_ICaL_Recovery, res.IV_ICaL_Recovery_Baseline, '-k', 'LineWidth', 2, 'MarkerSize', 8);
        subplot(2,6,[7,8]); plot(res.Exp_ICaL_Recovery_Baseline(:,1), res.Exp_ICaL_Recovery_Baseline(:,2), 'ok', res.t_ICaL_Recovery, res.IV_ICaL_Recovery_Baseline, '-k', res.t_ICaL_Recovery, res.IV_ICaL_Recovery_ISO, '-r', res.t_ICaL_Recovery, res.IV_ICaL_Recovery_CaMKII, '-b', 'LineWidth', 2, 'MarkerSize', 8);
        subplot(2,6,[7,8]); plot(res.t_ICaL_Recovery, res.IV_ICaL_Recovery_Baseline, '-k', res.t_ICaL_Recovery, res.IV_ICaL_Recovery_ISO, '-r', res.t_ICaL_Recovery, res.IV_ICaL_Recovery_CaMKII, '-b', 'LineWidth', 2, 'MarkerSize', 8);
        
        foptions = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [-Inf; 0], 'StartPoint', [-1.0; 50]);
        ftype = fittype('a * exp(-x / b) + 1', 'options', foptions);
        [fitR,gofR,outR] = fit(res.t_ICaL_Recovery,res.IV_ICaL_Recovery_Baseline,ftype);        
        mdl_tau_ICaL_Recovery_Baseline = fitR.b;        
        [fitR,gofR,outR] = fit(res.t_ICaL_Recovery,res.IV_ICaL_Recovery_CaMKII,ftype);        
        mdl_tau_ICaL_Recovery_CaMKII = fitR.b;        
        [fitR,gofR,outR] = fit(res.t_ICaL_Recovery,res.IV_ICaL_Recovery_ISO,ftype);        
        mdl_tau_ICaL_Recovery_ISO = fitR.b;
        
        subplot(2,6,9); bar([mdl_tau_ICaL_Recovery_Baseline, 77; mdl_tau_ICaL_Recovery_ISO, NaN; mdl_tau_ICaL_Recovery_CaMKII, NaN]);
        subplot(2,6,[10,11]); plot(Trace_Show_Baseline(:,1), Trace_Show_Baseline(:,2), '-k', Trace_Show_ISO(:,1), Trace_Show_ISO(:,2), '-r', Trace_Show_CaMKII(:,1), Trace_Show_CaMKII(:,2), '-b', 'LineWidth', 2);
        subplot(2,6,12); bar([mdl_tau_ICaL_ISO / mdl_tau_ICaL_Baseline, NaN; mdl_tau_ICaL_CaMKII / mdl_tau_ICaL_Baseline, 1.3]);
        
        resICaL = res;
        fid = fopen([outputfolder, 'Electrophysiology_IV_ICaL_ACT.txt'], 'w');
        fprintf(fid, 'variables=V \t Baseline \t PEND_Baseline \t ISO \t PEND_ISO \t CaMKII \t PEND_CaMKII\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f\n', [resICaL.V_ICaL_ACT, resICaL.IV_ICaL_ACT_Baseline, resICaL.IV_ICaL_PEND_Baseline, resICaL.IV_ICaL_ACT_ISO, resICaL.IV_ICaL_PEND_ISO, resICaL.IV_ICaL_ACT_CaMKII, resICaL.IV_ICaL_PEND_CaMKII]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_Exp_ICaL.txt'], 'w');
        fprintf(fid, 'variables=V \t Baseline\n');
        fprintf(fid, '%8.6f \t %8.6f\n', [resICaL.Exp_ICaL_IV_Baseline(:,1), resICaL.Exp_ICaL_IV_Baseline(:,2)]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_IV_ICaL_INACT.txt'], 'w');
        fprintf(fid, 'variables=V \t Baseline \t ISO \t CaMKII\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f\n', [resICaL.V_ICaL_INACT, resICaL.IV_ICaL_INACT_Baseline, resICaL.IV_ICaL_INACT_ISO, resICaL.IV_ICaL_INACT_CaMKII]');
        fclose(fid);        
    
        fid = fopen([outputfolder, 'Electrophysiology_Exp_ICaL_Inact.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_INACT_Baseline\n');
        fprintf(fid, '%8.6f \t %8.6f\n', [resICaL.Exp_ICaL_Inact_Baseline(:,1), resICaL.Exp_ICaL_Inact_Baseline(:,2)]');
        fclose(fid);
   
        fid = fopen([outputfolder, 'Electrophysiology_Comb_ICaL_Amp.txt'], 'w');
        fprintf(fid, 'variables=index \t Mdl \t Exp \t Std\n');
        fprintf(fid, '%f \t %8.6f \t %8.6f\n', [1, mv_iso ./ mv_base, 0.000, 0.000; 2, 0.000, 4.1, 0.31; 4, mv_camk ./ mv_base, 0.000, 0.000; 5, 0.000, 1.31, 0.09]');
        fclose(fid);
    
        fid = fopen([outputfolder, 'Electrophysiology_Comb_ICaL_Vshift.txt'], 'w');
        fprintf(fid, 'variables=index \t Mdl \t Exp\n');
        fprintf(fid, '%f \t %8.6f \t %8.6f\n', [1, abs(resICaL.V_ICaL_ACT(ind_iso) - resICaL.V_ICaL_ACT(ind_base)), 0.000; 2.2, 0.000, 15; 3.8, abs(mdl_vhalf_iso - mdl_vhalf_base), 0.000; 5, 0.000, 8.0]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_IV_ICaL_RECOVERY.txt'], 'w');
        fprintf(fid, 'variables=V \t Baseline \t ISO \t CaMKII\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f\n', [res.t_ICaL_Recovery, resICaL.IV_ICaL_Recovery_Baseline, resICaL.IV_ICaL_Recovery_ISO, resICaL.IV_ICaL_Recovery_CaMKII]');
        fclose(fid);
        
        fid = fopen([outputfolder, 'Electrophysiology_Comb_ICaL_RecTau.txt'], 'w');
        fprintf(fid, 'variables=index \t Mdl \t Exp \t Std\n');
        fprintf(fid, '%f \t %8.6f \t %8.6f \t %8.6f\n', [1, mdl_tau_ICaL_Recovery_Baseline, 0.000, 0.000; 2, 0.000, 69.7, 6.8; 4, mdl_tau_ICaL_Recovery_ISO, 0.000, 0.000; 7, mdl_tau_ICaL_Recovery_CaMKII, 0.000, 0.000]');
        fclose(fid);
        
        fid = fopen([outputfolder, 'Electrophysiology_Traces_ICaL.txt'], 'w');
        fprintf(fid, 'variables=t \t Baseline \t ISO \t CaMKII\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f\n', [(-10:0.5:200)', interp1(Trace_Show_Baseline(:,1), Trace_Show_Baseline(:,2), (-10:0.5:200)', 'spline'), interp1(Trace_Show_ISO(:,1), Trace_Show_ISO(:,2), (-10:0.5:200)', 'spline'), interp1(Trace_Show_CaMKII(:,1), Trace_Show_CaMKII(:,2), (-10:0.5:200)', 'spline')]');
        fclose(fid);
        
        fid = fopen([outputfolder, 'Electrophysiology_Comb_ICaL_InactTau.txt'], 'w');
        fprintf(fid, 'variables=index \t Mdl \t Exp \t Std\n');
        fprintf(fid, '%f \t %8.6f \t %8.6f \t %8.6f\n', [1, mdl_tau_ICaL_ISO ./ mdl_tau_ICaL_Baseline, 0.000, 0.000; 4, mdl_tau_ICaL_CaMKII ./ mdl_tau_ICaL_Baseline, 0.000, 0.000; 5, 0.000, 1.3, 0.000]');
        fclose(fid);
        
        pause;
        
    %% IKs
    % ===================================
        settings = [];
        settings.runSignalingPathway = 0;
        settings.ISO = 0;
        settings.K_o = 4.0;

        settings.INaB = 1;
        settings.ICaLB = 1;
        settings.INaCaB = 1;
        settings.INaKB = 1;
        settings.ITo2B = 1;
        settings.ICabB = 1;
        settings.applyVoltageClamp = 1;

        settings.fIKsP = 0.0;
        settings.fPLBP = 0.0;
        settings.fRyRP = 0.0;
        settings.fICaLP = 0.0;
        settings.fINaP = 0.0;   
        settings.fINaKP = 0.0;
        settings.fTnIP = 0.0;
        settings.fIKurP = 0;
        settings.k_JSR_Leak = 9E-4;
        
        vcp = [-50, 3000; 80, 6000; -25, 6500]; 
        [Traces, IV_IKs_Baseline_3000] = ComputeIVTraces(settings, vcp, (-20:10:70)', 2, 'iks', (2900:0.25:6400)', [6000.2, 6200]);

        vcp = [-50, 3000; 80, 3300; -25, 3800]; 
        [Traces, IV_IKs_Baseline_300] = ComputeIVTraces(settings, vcp, (-20:10:70)', 2, 'iks', (2900:0.25:3700)', [3300.2, 3500]);

        settings.fIKsP = 1.0;
        settings.fPLBP = 1.0;
        settings.fRyRP = 1.0;
        settings.fICaLP = 1.0;
        settings.fINaP = 1.0;   
        settings.fINaKP = 1.0;
        settings.fTnIP = 1.0;
        settings.fIKurP = 1.0;
        settings.k_JSR_Leak = 3.1E-3;
        
        vcp = [-50, 3000; 80, 6000; -25, 6500]; 
        [Traces, IV_IKs_ISO_3000] = ComputeIVTraces(settings, vcp, (-20:10:70)', 2, 'iks', (2900:0.25:6400)', [6000.2, 6200]);

        vcp = [-50, 3000; 80, 3300; -25, 3800]; 
        [Traces, IV_IKs_ISO_300] = ComputeIVTraces(settings, vcp, (-20:10:70)', 2, 'iks', (2900:0.25:3700)', [3300.2, 3500]);

        res.V_IKs = (-20:10:70)';
        res.IV_IKs_Baseline_3000 = IV_IKs_Baseline_3000'; % ./ max(IV_IKs_Baseline_3000);
        res.IV_IKs_Baseline_300 = IV_IKs_Baseline_300'; % ./ max(IV_IKs_Baseline_3000);
        res.IV_IKs_ISO_3000 = IV_IKs_ISO_3000';
        res.IV_IKs_ISO_300 = IV_IKs_ISO_300';

        res.Exp_IKs_Baseline_3000 = 2 * 2.4 * [0; 0.0061713; 0.037112; 0.11341; 0.24744; 0.433; 0.60825; 0.76291; 0.90722; 1];          % x 2 because EPI CELL!!
        res.Exp_IKs_Baseline_300 = 2 * 2.4 * [0; 0; 0; 0.012089; 0.024136; 0.060402; 0.12076; 0.16908; 0.21739; 0.25361];
        res.Exp_IKs_ISO_3000 = 2 * 7.2 * [0.040819; 0.098638; 0.21429; 0.36394; 0.54421; 0.70067; 0.82653; 0.91836;0.95578; 1];
        res.Exp_IKs_ISO_300 = 2 * 7.2 * [0.023903; 0.039843; 0.061765; 0.093631; 0.14344; 0.21118; 0.26696; 0.30481; 0.31676; 0.32274];
        
        figure;
            axh(1) = subplot(1,2,1);
                plot(res.V_IKs, res.Exp_IKs_Baseline_300, 'ok', res.V_IKs, res.Exp_IKs_ISO_300, 'sr', res.V_IKs, res.IV_IKs_Baseline_300, '-k', res.V_IKs, res.IV_IKs_ISO_300, '-r', 'LineWidth', 2, 'MarkerSize', 8);
            axh(2) = subplot(1,2,2);
                plot(res.V_IKs, res.Exp_IKs_Baseline_3000, 'ok', res.V_IKs, res.Exp_IKs_ISO_3000, 'sr', res.V_IKs, res.IV_IKs_Baseline_3000, '-k', res.V_IKs, res.IV_IKs_ISO_3000, '-r', 'LineWidth', 2, 'MarkerSize', 8);
            linkaxes(axh, 'y');
            ylim([-1 14]);
            
        resIKs = res;
        outputfolder = '';
        fid = fopen([outputfolder, 'Electrophysiology_IV_IKs_IV.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_300_Baseline \t EXP_300_Baseline \t IV_3000_Baseline \t EXP_3000_Baseline \t IV_300_ISO \t EXP_300_ISO \t IV_3000_ISO \t EXP_3000_ISO\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f\n', [resIKs.V_IKs, resIKs.IV_IKs_Baseline_300, resIKs.Exp_IKs_Baseline_300, resIKs.IV_IKs_Baseline_3000, resIKs.Exp_IKs_Baseline_3000, resIKs.IV_IKs_ISO_300, resIKs.Exp_IKs_ISO_300, resIKs.IV_IKs_ISO_3000, resIKs.Exp_IKs_ISO_3000]');
        fclose(fid);
            
    %% INa
    % =============================
        outputfolder = '';
        settings = [];
        settings.runSignalingPathway = 0;
        settings.ISO = 0.0;        
        settings.Na_o = 10;
        
        settings.ICaLB = 1.0;
        settings.IKsB = 1.0;    
        settings.IKrB = 1.0;
        settings.INaKB =1.0;
        settings.INaCaB = 0.0;
        settings.IKpB = 1.0;
        settings.IK1B = 1.0;
        settings.INabB = 0.0;
        settings.ITo1B = 1.0;
        settings.ITo2B = 0.0;
        settings.INaB = 0.0;
        settings.INaLB = 0.0;
        settings.IClB = 0.0;
        settings.IpCaB = 0.0;
        settings.CTKClB = 1.0;
        settings.CTNaClB = 0.0;
        settings.ICabB = 0.0;
        settings.IrelB = 0.0;
        settings.IupB = 0.0;
        settings.ItrB = 0.0;
        settings.IleakB = 0.0;
        settings.IdiffB = 0.0;
        
        settings.fICaLP = 0;
        settings.fRyRP = 0;
        settings.fIKsP = 0;
        settings.fPLBP = 0;
        settings.fINaP = 0;
        settings.fINaKP = 0;
        settings.fTnIP = 0;
        settings.fIKurP = 0;
        
        vcp_ACT = [-120, 1000; 10, 1050; -120, 1200];
        val_ACT = (-70:5:5)';
        ind_ACT = 2;
        [trace_ACT_NP, IV_ACT_NP] = ComputeIVTraces(settings, vcp_ACT, val_ACT, ind_ACT, 'ina', (950:0.1:1100)', [1000, 1050]);

        vcp_INACT = [-140, 1000; -50, 1500; -40, 1520; -140, 1820];
        val_INACT = [-120; -110; -100; -95; -90; -85; -80; -75; -70; -65; -60; -50];
        ind_INACT = 2;
        [trace_INACT_NP, IV_INACT_NP] = ComputeIVTraces(settings, vcp_INACT, val_INACT, ind_INACT, 'ina', (990:0.1:1550)', [1500, 1520]);
    
        settings.fICaLP = 1.0;
        settings.fRyRP = 1.0;
        settings.fIKsP = 1.0;
        settings.fPLBP = 1.0;
        settings.fINaP = 1.0;
        settings.fINaKP = 1.0;
        settings.fTnIP = 1.0;
        settings.fIKurP = 1.0;
        
        [trace_ACT_P, IV_ACT_P] = ComputeIVTraces(settings, vcp_ACT, val_ACT, ind_ACT, 'ina', (950:0.1:1100)', [1000, 1050]);
        [trace_INACT_P, IV_INACT_P] = ComputeIVTraces(settings, vcp_INACT, val_INACT, ind_INACT, 'ina', (990:0.1:1550)', [1500, 1520]);
    
        res.V_INa_ACT = val_ACT;
        res.IV_INa_ACT_Baseline = IV_ACT_NP';
        res.IV_INa_ACT_ISO = IV_ACT_P';
        res.V_INa_INACT = val_INACT;
        res.IV_INa_INACT_Baseline = IV_INACT_NP';
        res.IV_INa_INACT_ISO = IV_INACT_P';
        res.IV_GNa_ACT_Baseline = IV_ACT_NP'; % NOT CORRECT
        res.IV_GNa_ACT_ISO = IV_ACT_P';
        
        res.Exp_INa_NP_IV = [-0.01979798; -0.026262626; -0.029494949; -0.045656566; -0.090909091; -0.178181818; -0.346262626; -0.582222222; -0.811717172; -0.973333333; -1.002424242; -0.924848485; -0.756767677; -0.527272727; -0.255757576; 0.028686869];
        res.Exp_INa_P_IV = [-0.026262626; -0.032727273; -0.039191919; -0.065050505; -0.12969697; -0.294545455; -0.582222222; -0.928080808; -1.18989899; -1.325656566; -1.325656566; -1.228686869; -1.050909091; -0.818181818; -0.569292929; -0.288080808];
        res.Exp_INa_NP_Inact = [0.984962406; 0.992481203; 0.984962406; 0.962406015; 0.917293233; 0.827067669; 0.669172932; 0.481203008; 0.293233083; 0.157894737; 0.067669173; 0.015037594];
        res.Exp_INa_P_Inact = [1; 0.992481203; 0.954887218; 0.909774436; 0.827067669; 0.676691729; 0.481203008; 0.293233083; 0.135338346; 0.060150376; 0.022556391; 0.015037594];
        res.Exp_INa_NP_Act =  [0.005714286; 0.011428571; 0.011428571; 0.017142857; 0.028571429; 0.068571429; 0.148571429; 0.257142857; 0.411428571; 0.56; 0.691428571; 0.805714286; 0.891428571; 0.954285714; 0.982857143; 1];
        res.Exp_INa_P_Act = [0.005714286; 0.011428571; 0.011428571; 0.028571429; 0.062857143; 0.102857143; 0.211428571; 0.365714286; 0.52; 0.662857143; 0.777142857; 0.862857143; 0.92; 0.96; 0.988571429; 1]; 
        
        figure;
            subplot(1,2,1);
                plot(res.V_INa_ACT, res.IV_INa_ACT_Baseline ./ max(abs(res.IV_INa_ACT_Baseline)), '-k', res.V_INa_ACT, res.IV_INa_ACT_ISO ./ max(abs(res.IV_INa_ACT_Baseline)), '-r', res.V_INa_ACT, res.Exp_INa_NP_IV, 'ok', res.V_INa_ACT, res.Exp_INa_P_IV, 'sr', 'LineWidth', 2, 'MarkerSize', 8);
                
            subplot(1,2,2);
                plot(res.V_INa_INACT, res.IV_INa_INACT_Baseline ./ res.IV_INa_INACT_Baseline(1), '-k', res.V_INa_INACT, res.IV_INa_INACT_ISO ./ res.IV_INa_INACT_ISO(1), '-r', res.V_INa_INACT, res.Exp_INa_NP_Inact, 'ok', res.V_INa_INACT, res.Exp_INa_P_Inact, 'sr', 'LineWidth', 2, 'MarkerSize', 8);
        
        resINa = res;
        fid = fopen([outputfolder, 'Electrophysiology_IV_INa_ACT.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_ACT_Baseline \t IV_ACT_Baseline_GNa \t IV_ACT_ISO \t IV_ACT_ISO_GNa\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f\n', [resINa.V_INa_ACT, resINa.IV_INa_ACT_Baseline ./ max(abs(resINa.IV_INa_ACT_Baseline)), resINa.IV_GNa_ACT_Baseline ./ resINa.IV_GNa_ACT_Baseline(end), resINa.IV_INa_ACT_ISO ./ max(abs(resINa.IV_INa_ACT_Baseline)), resINa.IV_GNa_ACT_ISO ./ resINa.IV_GNa_ACT_ISO(end)]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_Exp_INa_ACT.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_ACT_Baseline \t IV_ACT_Baseline_GNa \t IV_ACT_ISO \t IV_ACT_ISO_GNa\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f \t %8.6f \t %8.6f\n', [resINa.V_INa_ACT, resINa.Exp_INa_NP_IV, resINa.Exp_INa_NP_Act, resINa.Exp_INa_P_IV, resINa.Exp_INa_P_Act]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_IV_INa_INACT.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_INACT_Baseline \t IV_INACT_ISO\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f\n', [resINa.V_INa_INACT, resINa.IV_INa_INACT_Baseline ./ resINa.IV_INa_INACT_Baseline(1), resINa.IV_INa_INACT_ISO ./ resINa.IV_INa_INACT_ISO(1)]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_Exp_INa_Inact.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_INACT_Baseline \t IV_INACT_ISO\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f\n', [resINa.V_INa_INACT, resINa.Exp_INa_NP_Inact, resINa.Exp_INa_P_Inact]');
        fclose(fid);

    %% RyR
    % ===============================================
% % %         [Tot_Rel_Baseline, Tot_Rel_ISO, Max_dRel_Baseline, Max_dRel_ISO, Max_dCa_Baseline, Max_dCa_ISO, res] = RyRPTest([1.9935; 0.5357]);
% % %         fid = fopen([outputfolder, 'Electrophysiology_RyR_Properties.txt'], 'w');
% % %         fprintf(fid, 'variables=index \t Mdl \t Exp\n');
% % %         fprintf(fid, '%d \t %8.6f \t %8.6f\n', [1, mean(Tot_Rel_ISO ./ Tot_Rel_Baseline), 0.000; 2, 0.000, 1.0; 4, mean(Max_dRel_ISO ./ Max_dRel_Baseline), 0.000; 5, 0.000, 1.867]');
% % %         fclose(fid);
        
    %% IKur
    % ===============================================
        settings = [];
        settings.runSignalingPathway = 0;
        settings.ISO = 0.0;                
        settings.K_o = 5;
        
        settings.ICaLB = 1.0;
        
        settings.fICaLP = 0;
        settings.fRyRP = 0;
        settings.fIKsP = 0;
        settings.fPLBP = 0;
        settings.fINaP = 0;
        settings.fINaKP = 0;
        settings.fTnIP = 0;
        settings.fIKurP = 0;        
        
        vcp_ACT = [-40, 1000; 30, 1080; -20, 1380; -40, 1500];
        val_ACT = (-30:10:60)';
        ind_ACT = 3;
        [trace_NP, IV_NP] = ComputeIVTraces(settings, vcp_ACT, val_ACT, ind_ACT, 'ikp', (990:0.5:1400)', [1100, 1379]);

        settings.fICaLP = 1.0;
        settings.fRyRP = 1.0;
        settings.fIKsP = 1.0;
        settings.fPLBP = 1.0;
        settings.fINaP = 1.0;
        settings.fINaKP = 1.0;
        settings.fTnIP = 1.0;
        settings.fIKurP = 1.0;
        
        [trace_P, IV_P] = ComputeIVTraces(settings, vcp_ACT, val_ACT, ind_ACT, 'ikp', (990:0.5:1400)', [1100, 1395]);
        
        res.V_IKp = val_ACT;
        res.IV_IKp_Baseline = IV_NP';
        res.IV_IKp_ISO = IV_P';
        
        res.Exp_IKp_Baseline = [(-20:10:50)',[0.005960265; 0.068543046; 0.124172185; 0.145033113; 0.179801325; 0.263245033; 0.381456954; 0.492715232]];
%         res.Exp_IKp_ISO = [(-20:10:50)',[0.054635762; 0.158940397; 0.145033113; 0.360596026; 0.57615894; 0.882119205; 1.18807947; 1.438410596]];
        res.Exp_IKp_ISO = [(-20:10:50)',[0.027475248; 0.079455446; 0.141831683; 0.245792079; 0.412128713; 0.640841584; 0.864356436; 1.186633663]];

        figure;
            plot(res.Exp_IKp_Baseline(:,1), res.Exp_IKp_Baseline(:,2), 'ob', res.Exp_IKp_ISO(:,1), res.Exp_IKp_ISO(:,2), 'sr', res.V_IKp, res.IV_IKp_Baseline, '-b', res.V_IKp, res.IV_IKp_ISO, '-r', 'LineWidth', 2);        
        
        resIKur = res;
        fid = fopen([outputfolder, 'Electrophysiology_IV_IKur.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_ACT_Baseline \t IV_ACT_ISO\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f\n', [resIKur.V_IKp, resIKur.IV_IKp_Baseline, resIKur.IV_IKp_ISO]');
        fclose(fid);

        fid = fopen([outputfolder, 'Electrophysiology_Exp_IKur.txt'], 'w');
        fprintf(fid, 'variables=V \t IV_ACT_Baseline \t IV_ACT_ISO\n');
        fprintf(fid, '%8.6f \t %8.6f \t %8.6f\n', [resIKur.Exp_IKp_Baseline(:,:), resIKur.Exp_IKp_ISO(:,2)]');
        fclose(fid);
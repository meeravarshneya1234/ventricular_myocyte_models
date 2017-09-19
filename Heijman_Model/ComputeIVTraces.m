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
% *   When called with a settings structure, voltage-clamp protocol and   *
% *   voltage-clamp values (as well as the index in the vcp that should   *
% *   be changed) and the name of the current of interest, this function  *
% *   returns the I-V relationship of that current at a certain time of   *
% *   interest. In addition a set of traces is returned (using the        *
% *   specified interpolation)                                            *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   mainHRdBA.m,                                                        *
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

function [Traces, Peaks, X0, State, Ti, currents] = ComputeIVTraces(settings, vcp, vvalues, vcpindex, current, interp, toi)
    nu = 1;
    
    Traces = zeros(length(interp), nu*length(vvalues));
    Peaks = zeros(1, length(vvalues));
    
    h4 = waitbar(0,' Running Voltage-Clamp Protocol ...');
    for i=1:length(vvalues)
        vcpcopy = vcp;
        toicopy = toi;
        
        if vcpindex > 0
            vcpcopy(vcpindex,1) = vvalues(i);
        else
            a = -vcpindex(1);
            vcpcopy(a:end, 2) = vcpcopy(a:end,2) + vvalues(i) * ones(size(vcpcopy,1)-a+1,1);
            toicopy(1) = toicopy(1) + vvalues(i);
            toicopy(2) = toicopy(2) + vvalues(i);
        end
        
        settings.vcp = vcpcopy;
        settings.applyVoltageClamp = 1;
        
        [currents, State, Ti] = mainHRdBA(settings);     
        X0 = State(end,:);
        
        yvals = [];
        if isa(current, 'char')
            if isfield(currents, current), yvals = getfield(currents, current); end
        else
            yvals = State(:, current);
        end
        if strcmpi(current, 'RyR')
            nu = 4;
            yvals = [currents.ical, State(:, 2), State(:, 5), currents.irel]; 
        end
        
        [Ti, m] = unique(Ti);
        yvals = yvals(m,:);
        
        [v, si] = min(abs(Ti - toicopy(1)));
        [v, ei] = min(abs(Ti - toicopy(2)));
        
        ind = getPeaks([Ti(si+1:ei-1), abs(yvals(si+1:ei-1,1))], [10, NaN, 1]);
        if length(ind) ~= 1
            %disp('Multiple or no peak(s) found!');
            [v, ind] = max(abs(yvals(si+1:ei-1,1)));
            if ~isempty(ind)
                Peaks(i) = yvals(si+1+ind,1);
            else
                Peaks(i) = NaN;
            end
        else
            Peaks(i) = yvals(si+1+ind,1);
        end
        
        yinterp = interp1(Ti, yvals, interp);
        
        if length(vvalues) > 1            
            Traces(:, nu*(i-1)+1:nu*i) = yinterp;
        else
            Traces = [Ti, yvals];
        end
        
        waitbar(i/length(vvalues),h4)
    end
    
    close(h4);
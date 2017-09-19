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
% *   When called with an Nx2 matrix and some options, this function will *
% *   detect the peaks in the signal (first column time axis, second      *
% *   column data. The first option specifies the number of points to the *
% *   left and right that should be smaller than the current point to     *
% *   count as a peak. The second options specifies the minimum amplitude *
% *   (as a percentage of maximum amplitude in signal). The final option  *
% *   specifies the minimum distance between peaks.                       *
% *                                                                       *
% *   This m-file depends on the following files:                         *
% *   -                                                                   *
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

function indices = getPeaks(data, options)
    indices = [];
    if nargin == 1
        options = [1, 0.5, 1]; %range, max val perc, dist
    end
    
    maxval = max(data(:,2));
    
    for i=options(1)+1:size(data,1)-options(1)
        if isnan(options(2)) || (~isnan(options(2)) && data(i,2) > options(2) * maxval)
            if  data(i,2) >= max(data(i-options(1):i-1,2))
                if data(i,2) >= max(data(i+1:i+options(1),2))
                    if length(indices)==0
                        indices = [indices, i];
                    elseif (data(i,1) - data(indices(end),1)) >= (options(3) / 1000)
                        indices = [indices, i];
                    end
                end
            end
        end
    end
                        
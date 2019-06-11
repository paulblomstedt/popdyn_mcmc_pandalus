function phiR = growth_recruits2(Gk, Linfsigma, Linf, Gtime, t0, midpoints, I)
% Population state transition due to growth as described in Samu's GPDM manuscript
% phiG = expected proportions of individuals in each length class after growth
% nG = no. of individuals in each length group after growth

% Copyright (c) 2013 Jarno Vanhatalo
% Copyright (c) 2014 Paul Blomstedt

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


    LC = length(midpoints);              % number of length intervals
    
    % Growth matrix: normally distributed growths from each length class.
    % I: end points of length classes
    % Von Bertalanffy Growth model: expected length from size class k
    ELR = Linf*(1-exp(-Gk*(Gtime+t0)));
    SLR = sqrt(Linfsigma*Linfsigma*(1-exp(-Gk*(Gtime+t0)*2)));
    tmp = normcdf((I-ELR)/SLR);
    phiR = ( tmp(2:end)-tmp(1:end-1) )./( tmp(end)-tmp(1) );        % Length distribution of recruits, defined using the growth equation


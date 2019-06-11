function nG = growth_expected_num2(popstate, Gk, Linfsigma, Linf, Gtime, midpoints, I)
% Population state transition due to growth as described in Samu's GPDM manuscript
% phiG = expected proportions of individuals in each length class after growth
% nG = no. of individuals in each length group after growth

% Copyright (c) 2014 Paul Blomstedt

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


    LC = length(midpoints);              % number of length intervals
    
    % Growth matrix: normally distibuted growths from each length class.
    % I: end points of length classes
    % Von Bertalanffy Growth model: expected length from size class k
    Eg = (Linf-midpoints)*(1-exp(-Gk*Gtime))+midpoints;             % In time, the length distribution of a cohort converges; 
    Lsigma = sqrt(Linfsigma*Linfsigma*(1-exp(-Gk*Gtime*2)));        % Standard deviation of growths over a  time step, Linfsigma: standard deviation of lengths of old fish  
    g_ij = zeros(LC);
    for i=1:LC        
        tmp = normcdf((I-Eg(i))/Lsigma);
        g_ij(i,:) = ( tmp(2:end)-tmp(1:end-1) )./( tmp(end)-tmp(1) );    % probability of moving from class i to j
    end
    nG = (popstate'*g_ij)';                                                 % state of population after growth
    
    %NB! By the law of large numbers, the variance of the distribution of 
    % individuals in each class, relative to the total numbers, will become 
    % negligible, and therefore expected numbers will be used instead of 
    % random draws from a multinomial distribution.

    
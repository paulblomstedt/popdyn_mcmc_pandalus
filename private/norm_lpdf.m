function lp = norm_lpdf(x,mu,s2)
% Gaussian log probability density (LPDF)

% norm_lpdf(x,mu,s2) calculates the log probability density of Gaussian
% distribution at x with mean mu and variance s2

% Copyright (c) 2014 Paul Blomstedt

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

lp = -0.5*(log(2*pi) +log(s2) +1./s2.*(x-mu).^2);


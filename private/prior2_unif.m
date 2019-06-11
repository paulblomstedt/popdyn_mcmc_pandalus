function p = prior2_unif(varargin)
%PRIOR_UNIF  Uniform prior structure     
%       
%  Description
%    P = PRIOR_UNIF('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) creates
%    uniform prior structure. in which the named parameters have the
%    specified values. Any unspecified parameters are set to default
%    values. 
%    
%  See also
%    PRIOR_*

% Copyright (c) 2009, 2013 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_UNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('minp',-inf, @(x) isscalar(x) && isreal(x));
  ip.addParamValue('minp_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('maxp',inf, @(x) isscalar(x) && isreal(x));
  ip.addParamValue('maxp_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  % Initialize parameters
  if init || ~ismember('minp',ip.UsingDefaults)
      p.minp = ip.Results.minp;
  end
  if init || ~ismember('maxp',ip.UsingDefaults)
      p.maxp = ip.Results.maxp;
  end

  % Initialize prior structure
  if init
      p.p=[];
  end
  if init || ~ismember('minp_prior',ip.UsingDefaults)
      p.p.mu=ip.Results.minp_prior;
  end
  if init || ~ismember('maxp_prior',ip.UsingDefaults)
      p.p.s2=ip.Results.maxp_prior;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_unif_pak;
    p.fh.unpak = @prior_unif_unpak;
    p.fh.lp = @prior_unif_lp;
    p.fh.lpg = @prior_unif_lpg;
    p.fh.recappend = @prior_unif_recappend;
    p.fh.rnd = @prior_unif_rnd;
  end
  
end

function [w, s] = prior_unif_pak(p, w)
  w=[];
  s={};
  
  if ~isempty(p.p.minp)
    w = p.minp;
    s=[s; 'Uniform.minp'];
  end
  if ~isempty(p.p.maxp)
    w = [w p.maxp];
    s=[s; 'Uniform.maxp'];
  end

end

function [p, w] = prior_unif_unpak(p, w)
  if ~isempty(p.p.minp)
      i1=1;
      p.minp = w(i1);
      w = w(i1+1:end);
  end
  if ~isempty(p.p.maxp)
      i1=1;
      p.maxp = w(i1);
      w = w(i1+1:end);
  end
end

function lp = prior_unif_lp(x, p)
  
  ind = x<p.maxp & x>p.minp;
  lp = -inf*ones(size(ind));
  lp(ind) = 0;
%   if x<p.maxp && x>p.minp
%       lp = 0;
%   else
%       lp=-inf;
%   end
end

function lpg = prior_unif_lpg(x, p)
  lpg = zeros(size(x));
end

function rec = prior_unif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  if ~isempty(p.p.minp)
    rec.minp(ri,:) = p.minp;
  end
  if ~isempty(p.p.maxp)
    rec.maxp(ri,:) = p.maxp;
  end
end

function x = prior_unif_rnd(p)
% draw sample from the prior
    
  x = p.minp + (p.maxp-p.minp).*rand;
end



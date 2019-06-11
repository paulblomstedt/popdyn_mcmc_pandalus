function p = prior2_beta(varargin)
%PRIOR_LOGGAUSSIAN  beta prior structure     
%       
%  Description
%    P = PRIOR_BETA('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates Beta prior structure in which the named parameters have the
%    specified values. Any unspecified parameters are set to default
%    values. 
%    
%    P = PRIOR_BETA(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a prior structure with the named parameters altered
%    with the specified values.
%
%    Parameters for Beta prior [default]
%      a       - [1]
%      b       - [1]
%      a_prior - prior for a [prior_fixed]
%      b_prior - prior for b [prior_fixed]
%
%  See also
%    PRIOR_*

% Copyright (c) 2000-2001,2010 Aki Vehtari
% Copyright (c) 2010 Jaakko Riihim�ki
% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_BETA';
  ip.addOptional('p', [], @isstruct);
  ip.addParamValue('a',1, @(x) isscalar(x));
  ip.addParamValue('a_prior',[], @(x) isstruct(x) || isempty(x));
  ip.addParamValue('b',1, @(x) isscalar(x) && x>0);
  ip.addParamValue('b_prior',[], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Beta';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Beta')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('a',ip.UsingDefaults)
    p.a = ip.Results.a;
  end
  if init || ~ismember('b',ip.UsingDefaults)
    p.b = ip.Results.b;
  end
  % Initialize prior structure
  if init
    p.p=[];
  end
  if init || ~ismember('a_prior',ip.UsingDefaults)
    p.p.a=ip.Results.a_prior;
  end
  if init || ~ismember('b_prior',ip.UsingDefaults)
    p.p.b=ip.Results.b_prior;
  end

  if init
    % set functions
    p.fh.pak = @prior_beta_pak;
    p.fh.unpak = @prior_beta_unpak;
    p.fh.lp = @prior_beta_lp;
    p.fh.lpg = @prior_beta_lpg;
    p.fh.recappend = @prior_beta_recappend;
    p.fh.rnd = @prior_beta_rnd;
  end

end

function [w, s] = prior_beta_pak(p)
  
  w=[];
  s={};
  if ~isempty(p.p.a)
    w = p.a;
    s=[s; 'beta.a'];
  end
  if ~isempty(p.p.b)
    w = [w p.b];
    s=[s; 'beta.b'];
  end
end

function [p, w] = prior_beta_unpak(p, w)

  if ~isempty(p.p.a)
    i1=1;
    p.a = w(i1);
    w = w(i1+1:end);
  end
  if ~isempty(p.p.b)
    i1=1;
    p.b = w(i1);
    w = w(i1+1:end);
  end
end

function lp = prior_beta_lp(x, p)
  
  if any(x<0 | x>1)
      lp = -inf;
  else
      lp= (p.a-1).*log(x) +(p.b-1).*log(1-x) - betaln(p.a,p.b);
  end

  if ~isempty(p.p.a)
    lp = lp + p.p.a.fh.lp(p.a, p.p.a);
  end
  if ~isempty(p.p.b)
    lp = lp + p.p.b.fh.lp(p.b, p.p.b);
  end
end

function lpg = prior_beta_lpg(x, p)
  
  lpg = (p.a-1)./x - (p.b-1)./(1-x);
  
  if ~isempty(p.p.a)
      error('gradient with respect to a is not implemented.')
    lpga = [];
    lpg = [lpg lpga];
  end
  if ~isempty(p.p.b)
      error('gradient with respect to b is not implemented.')
    lpgb = [];
    lpg = [lpg lpgb];
  end
end

function rec = prior_beta_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  if ~isempty(p.p.a)
    rec.a(ri,:) = p.a;
  end
  if ~isempty(p.p.b)
    rec.b(ri,:) = p.b;
  end
end    

function x = prior_beta_rnd(p)
% draw sample from the prior
    
  x = betarand(p.a,p.b);
end    




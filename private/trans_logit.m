function t = trans_logit(varargin)
% TRANS_LOGIT   logit transformation function
%
%  Description
%    P = TRANS_LOGIT('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates logit transformation structure in which the named parameters
%    have the specified values. Any unspecified parameters are set to
%    default values. 
%    
%    P = TRANS_LOGIT(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a transformation structure with the named parameters altered
%    with the specified values.
%
%    Parameters for truncated-Gaussian prior [default]
%      minp     - minimum of the non transformed parameter [0]
%      maxp     - minimum of the non transformed parameter [0]
%
%    the transformation is tx = logit( (x - t.minp)./(t.maxp-t.minp) );

% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'TRANS_LOGIT';
  ip.addOptional('t', [], @isstruct);
  ip.addParamValue('minp',0, @(x) isscalar(x) && isreal(x));
  ip.addParamValue('maxp',1, @(x) isscalar(x) && isreal(x));
  ip.parse(varargin{:});
  t=ip.Results.t;
  
  if isempty(t)
    init=true;
    t.type = 'trans-logit';
  else
    if ~isfield(t,'type') && ~isequal(t.type,'trans-logit')
      error('First argument does not seem to be a valid transformation structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('minp',ip.UsingDefaults)
    t.minp = ip.Results.minp;
  end
  if init || ~ismember('maxp',ip.UsingDefaults)
    t.maxp = ip.Results.maxp;
  end

  if init
    % set functions
    t.fh.trans = @trans_logit_trans;
    t.fh.invtrans = @trans_logit_invtrans;
    t.fh.jit = @trans_logit_jit;
    t.fh.ljit = @trans_logit_ljit;
  end
end

function tx = trans_logit_trans(x, t)
    tx = logit( (x - t.minp)./(t.maxp-t.minp) );
end

function x = trans_logit_invtrans(tx, t)
    x = t.minp + (t.maxp-t.minp)*logitinv(tx);
end

function j = trans_logit_jit(tx, t)
    % jacobian of the inverse transformation
    maxcut = -log(eps);
    mincut = -log(1/realmin - 1);
    tt = max(min(tx,maxcut),mincut);
    j = (t.maxp-t.minp)*exp(-tt)./(1+exp(-tt)).^2;
    
end

function j = trans_logit_ljit(tx, t)
    % log jacobian of the inverse transformation
    maxcut = -log(eps);
    mincut = -log(1/realmin - 1);
    tt = max(min(tx,maxcut),mincut);
    j = log(t.maxp-t.minp) - tt -2.*log(1+exp(-tt));    
end




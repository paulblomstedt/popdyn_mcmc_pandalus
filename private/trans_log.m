function t = trans_log(varargin)
% TRANS_LOG   log transformation function
%
%  Description
%    P = TRANS_LOG('PARAM1', VALUE1, 'PARAM2', VALUE2, ...) 
%    creates log transformation structure in which the named parameters
%    have the specified values. Any unspecified parameters are set to
%    default values. 
%    
%    P = TRANS_LOGIT(P, 'PARAM1', VALUE1, 'PARAM2', VALUE2, ...)
%    modify a transformation structure with the named parameters altered
%    with the specified values.
%
%    Parameters for truncated-Gaussian prior [default]
%      minp     - minimum of the non transformed parameter [0]
%
%    the transformation is tx = log( x - t.minp );

% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'TRANS_LOG';
  ip.addOptional('t', [], @isstruct);
  ip.addParamValue('minp',0, @(x) isscalar(x) && isreal(x));
  ip.parse(varargin{:});
  t=ip.Results.t;
  
  if isempty(t)
    init=true;
    t.type = 'trans-log';
  else
    if ~isfield(t,'type') && ~isequal(t.type,'trans-log')
      error('First argument does not seem to be a valid transformation structure')
    end
    init=false;
  end

  % Initialize parameters
  if init || ~ismember('minp',ip.UsingDefaults)
    t.minp = ip.Results.minp;
  end

  if init
    % set functions
    t.fh.trans = @trans_log_trans;
    t.fh.invtrans = @trans_log_invtrans;
    t.fh.jit = @trans_log_jit;
    t.fh.ljit = @trans_log_ljit;
  end
end

function tx = trans_log_trans(x, t)
    tx = log( x - t.minp );
end

function x = trans_log_invtrans(tx, t)
    x = t.minp + exp(tx);
end

function j = trans_log_jit(tx, ~)
    % jacobian of the inverse transformation
    j = exp(tx);
    
end

function j = trans_log_ljit(tx, ~)
    % log jacobian of the inverse transformation
    j = tx;
    
end

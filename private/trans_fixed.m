function t = trans_fixed(varargin)
% TRANS_FIXED   fixed transformation function
%
%  Description
%    P = TRANS_FIXED 
%    creates a transformation structure which returns the parameter itself.
%    The transformation is tx = x;

% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'TRANS_LOG';
  ip.addOptional('t', [], @isstruct);
  ip.parse(varargin{:});
  t=ip.Results.t;
  
  if isempty(t)
    init=true;
    t.type = 'trans-fixed';
  else
    if ~isfield(t,'type') && ~isequal(t.type,'trans-fixed')
      error('First argument does not seem to be a valid transformation structure')
    end
    init=false;
  end

  if init
    % set functions
    t.fh.trans = @trans_fixed_trans;
    t.fh.invtrans = @trans_fixed_invtrans;
    t.fh.jit = @trans_fixed_jit;
    t.fh.ljit = @trans_fixed_ljit;
  end
end

function tx = trans_fixed_trans(x, ~)
    tx = x;
end

function x = trans_fixed_invtrans(tx, ~)
    x = tx;
end

function j = trans_fixed_jit(~, ~)
    % jacobian of the inverse transformation
    j = 1;
    
end

function j = trans_fixed_ljit(~, ~)
    % jacobian of the inverse transformation
    j = 0;
    
end

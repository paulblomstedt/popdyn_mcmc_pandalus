function [state, params, rec] = popdyn_mc2_pandalus(model, opt, rec, varargin)
% This function implements a Markov chain sampler for the pandalus 
% population dynamics model

% Copyright (c) 2014 Paul Blomstedt
% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


%% Start the inference
% ====================
T = model.T;

if isempty(rec)
    % Sample initial parameters and states from the prior
    [state, params] = model.fh.initialize(1);
    for t=1:T+1
        state = model.fh.state_forward_sample(t, state, params, model, 'state_forward_sample');
    end
    rec_s = [];
    rec_p = [];
else
    rec_s = rec.s;
    rec_p = rec.p;
    state = rec.state;
    params = rec.params;
end
    
pw = model.fh.pack_params(params);
nparams = length(pw);
% nN = length(state.N(:,1));

% if options structure is empty set the default options
if isempty(opt)
    nsamp = 10000;
    repeat = 10;
    prop_std_s = 0.1; 
    prop_std = 0.01*eye(nparams);
else
    nsamp = opt.nsamp;
    repeat = opt.repeat;
    prop_std = opt.prop_std;
    prop_std_s = opt.prop_std_s;
end

if isfield(opt, 'state_repeat')
    state_repeat = opt.state_repeat;
else
    state_repeat = 1;
end
if isfield(opt, 'param_repeat')
    param_repeat = opt.param_repeat;
else
    param_repeat = 1;
end

acc_p=0;
acc_s=0;

%%
for i1 = 1:nsamp % loop over iterations
    
    for k1=1:repeat % loop for thinning

        % sample states
        % ------------------
        acc=0;
        for ii=1:state_repeat
            for t=1:T % loop over time steps
                % Evaluate log target density at current point
                lp_old = model.fh.state_logprior(t, state, params, model) + ...                             
                    model.fh.state_logprior(t+1, state, params, model) + ...
                    model.fh.log_likelihood(t, state, params, model, 'log_likelihood', varargin{:});
                state_new = state;

                % LFmax proposal
                state_new.LFmax(:,t) = state.LFmax(:,t) + prop_std_s*randn;                 % GENERATE proposal for LFmax as random walk

                % LR proposal
                state_new.LR(:,t) = state.LR(:,t) + prop_std_s*randn;                       % GENERATE proposal for LR as random walk

                % Evaluate log target density at proposed point
                lp_new = model.fh.state_logprior(t, state_new, params, model) + ...             
                    model.fh.state_logprior(t+1, state_new, params, model) + ...
                    model.fh.log_likelihood(t, state_new, params, model, 'log_likelihood', varargin{:});

                % M-H test
                if rand < exp(lp_new - lp_old)              
                    if t==1
                        state = state_new;
                    else
                        state = model.fh.update_state(t, state, state_new, params, model);  % function which updates the remaining state variables deterministically
                    end
                    acc = acc+1;
                    lp_old = lp_new;
                end

            end % end loop over time steps
        end
        acc_s = acc_s+acc/(T*state_repeat);
        
        % sample parameters
        % ------------------
        for jj=1:param_repeat
            % Draw a proposal
            pw_new = pw + prop_std*randn(nparams,1);
            % calculate densitites
            lp_old = -nlog_density(pw, params, state, model, T, varargin{:});
            lp_new = -nlog_density(pw_new, params, state, model, T, varargin{:});
            % Accept / reject
            if rand < exp(lp_new - lp_old)
                pw = pw_new;
                params = model.fh.unpack_params(pw, params);
                acc_p = acc_p+1./param_repeat;
            end
        end
        
    end % end loop for thinning
    
    % record
    rec_s(:,end+1) = sample_pack(state);
    rec_p(:,end+1) = pw;
        
    if opt.display > 0                      
        fprintf('%d  err: %.2f acc_s: %.2f  acc_p: %.2f \n', size(rec_s,2), lp_old, acc_s./i1/repeat, acc_p./i1/repeat)
    end
    
    if opt.display > 1 && mod(i1,100)==0   
        params_t = model.fh.unpack_params(rec_p,params);
        state_t = sample_pack(state,rec_s);
        
        % plot parameter trajectories
        fnames = fieldnames(params.p);
        spi = 0;
        for jjj = 1:ceil(length(fnames)./5)
            figure(jjj)
            for iii=1:5
                spi = spi+1;
                subplot(5,1,iii), plot(params_t.(fnames{spi})), title(fnames{spi},'Interpreter','none')
                if spi==length(fnames)
                    break
                end
            end
        end
        
        % plot total population trajectory
         figure(ceil(length(fnames)./5)+1)
        for j1=2:T 
            subplot(5,6,j1-1), plot(squeeze(sum(state_t.N(:,j1,:,:),1))), title(['Tot. pop. size, t=',num2str(j1-1)])
        end
    end
    
    if opt.save == 1 && mod(i1,100)==0          % saves on every 100th iteration
        rec.s = rec_s;                          % saves all iterations as stacked vectors   
        rec.p = rec_p;
        rec.acc = [acc_p acc_s]./i1/repeat;
        rec.params = params;
        rec.state = state;                      % saves the last iteration as a structure   
        
        if isfield(opt, 'savename')
            save([opt.savename '_mcmc'], 'rec', 'model');
        else
            save([model.name '_mcmc'], 'rec', 'model');
        end
    end
    
end % end loop over iterations

rec.s = rec_s;
rec.p = rec_p;
rec.params = params;
rec.state = state;
rec.acc = [acc_p acc_s]./i1/repeat;

params = model.fh.unpack_params(rec_p,params);
state = sample_pack(state,rec_s);

end

function lp = nlog_density(pw, params, state, model, T, varargin)
% The negative log density of parameters
%  ! This has to be negative so that, e.g., slice-sampling works !

pw = pw(:);   % force pw to be column vector
params = model.fh.unpack_params(pw, params);
% Calculate log density
lp = model.fh.param_logprior(params);
for t=1:T 
    lp = lp + model.fh.state_logprior(t, state, params, model, varargin{:}) + ...
         model.fh.log_likelihood(t, state, params, model, 'log_likelihood', varargin{:});
end
lp = -lp;

end



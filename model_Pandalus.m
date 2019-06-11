function model = model_Pandalus(varargin)

% Copyright (c) 2014 Paul Blomstedt 
% Copyright (c) 2013 Jarno Vanhatalo 

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% About model dimensions
%  State variables are stored in matrices for which dimensions are:
%    [LC, T, FL, S], where
%  LC = length classes
%  T  = timesteps (=years*seasons)
%  FL = fleets
%  S  = Monte Carlo samples
%
% Parameters are stored in matrices of size:
%   [Np S], where 
% Np = length of parameter vector
% S  = number of samples

% NAMING CONVENTION
%
% N:    state of the population, numbers of individuals in each length class at the end of a time step
% n:    number of individuals per length class in a specified subset of the population
%   nG:     number of individuals per length class after growth
%   nS:     number of survivors per length class at the end of a time step
%   nC:     number of caught individuals per length-class
%   nR:     number of recruits per length class at the end of a time step
%   NB! nS+nR = N
% phi:  relative proportion of individuals per length class in a specified subset of the population 
%   phiR:   (normalized) length distribution of recruits
%   NB! sum(phiR)=1
% pi:   expected proportions of the entire population
%   piS:     probability that an individual in a given length class survives
%   piC:     probability that an individual in a given length class dies of fishing
%   piD:     probability that an individual in a given length class dies of reasons other than fishing ("natural" causes)
%   NB! sum([piS(i) piC(i,t) piD(i)])=1, i = 1,...,LC;

ip=inputParser;
ip.FunctionName = 'model_Pandalus';
ip.addParamValue('lengthIntervals', [], @(x) isreal(x) && all(x>=0))
ip.addParamValue('T', [], @(x) isreal(x) && isscalar(x) && x>0)
ip.addParamValue('FL', 1, @(x) isreal(x) && isscalar(x) && x>0)

ip.parse(varargin{:});
I=ip.Results.lengthIntervals;
T=ip.Results.T;
FL=ip.Results.FL;

model.name = 'Pandalus';
model.I = I;           % length intervals
model.T = T;           % number of time steps
model.FL = FL;         % number of fleets

model.fh.initialize = @initialize;
model.fh.sample_params_from_prior = @sample_params_from_prior;
model.fh.pack_params = @pack_params;
model.fh.unpack_params = @unpack_params;
model.fh.take_nth_params = @take_nth_params;
model.fh.take_nth_state = @take_nth_state;
model.fh.state_forward_sample = @state_forward_sample;
model.fh.log_likelihood = @log_likelihood;
model.fh.param_logprior = @param_logprior;
model.fh.state_logprior = @state_logprior;
model.fh.update_state = @update_state;

end

function [state params] = initialize(nsamp)

% Initialize the state structure
% ===============================
% ===============================
state.N = [];       % state of population at the end of each time period 
state.LFmax = [];   % stochasticity term related to fishing mortality
state.LR = [];      % stochasticity term related to stock-recruitment
state.nC = [];      % number of caught individuals per length-class
state.nG = [];      % state of population after growth
state.piS = [];     % probability that an individual in a given length class survives
state.nR = [];      % length distribution of recruits; monitored for convenience (not needed for forward sampling)
state.Fmax = [];    % maximum fishing mortality; monitored for convenience (not needed for forward sampling)

% Initialize params structure
% ===========================
% ===========================

% Some fixed parameters passed on as hyperparameters to priors
minRautoc = 0.1; maxRautoc = 0.9;   % min and max autoregression coefficient for recruitment
minFautoc = 0.1; maxFautoc = 0.9;   % min and max autoregression coefficient for fishing mortality

% Growth, priors
% ======================================================================= %
% priors for growth parameters
params.p.logitLinf = prior2_gaussian;                           
params.t.logitLinf = trans_fixed;                           
params.p.logitLinfsigma = prior2_gaussian;                      
params.t.logitLinfsigma = trans_fixed;                      
params.p.logitGk = prior2_gaussian;       
params.t.logitGk = trans_fixed;                           

% Mortality and survival, priors
% ======================================================================= %
% Assumption: natural mortality exceeds fishing mortality (stated in the stock annex, p.5)

% Natural mortality, priors
% ----------------------------------------------------------------------- %
params.p.Mconstant =...                                                     % "mean" natural mortality 
    prior2_loggaussian('mu', log(0.75), 's2', 0.001);                       %  The instananeous natural mortality used in Nielsen et al. 2013 is 0.75
params.t.Mconstant = trans_log;                         

% Fishing mortality, priors
% ----------------------------------------------------------------------- %

% "mean" fishing mortality
params.p.Fconstant = prior2_loggaussian('mu', log(0.4), 's2', 0.05);        % over-all fishing mortality over all time steps 
params.t.Fconstant = trans_log;

% autocorrelation parameters
params.p.F_autoc = prior2_unif('minp',minFautoc, 'maxp', maxFautoc);        % autoregression coefficient
params.t.F_autoc = trans_logit('minp', minFautoc, 'maxp', maxFautoc);
params.p.CVF = prior2_beta('a',1,'b',1);                                        % coefficient of variation for fishing mortality...
params.t.CVF = trans_logit;                                                     %   the original prior was CVF = 0.001 + 0.999*rand(1,1,FL);    % dunif(0.001,1)

% Selection curve parameters
params.p.L50 = prior2_loggaussian('mu',log(18), 's2', 0.005);               % length at 50% selectivity. 
params.t.L50 = trans_log;
params.p.betaGsel = prior2_loggaussian('mu',log(0.3), 's2', 0.1);           % softness of the selection curve. 
params.t.betaGsel = trans_log;


% Reproduction, priors
% ======================================================================= %

% Egg production, priors
% ----------------------------------------------------------------------- %

% Parameters for fecundity, no. of eggs = af*length^bf
params.p.logaf = prior2_gaussian('mu', -2, 's2', 1);  % NB! the nontation \alpha_f = logaf is used in the manuscript
params.t.logaf = trans_fixed;
params.p.bf = prior2_gaussian('mu', 3, 's2', 0.4);      %                 \beta_f = bf
params.t.bf = trans_fixed;

% Recruitment, priors
% ----------------------------------------------------------------------- %
% alpha = The survival at the limit of zero population
% K = carrying capacity (max no. of recruits relative to the initial population)
params.p.logitalpha = prior2_gaussian('mu', -3, 's2', 1);   
params.t.logitalpha = trans_fixed;
% The hyperparameters for K were obtained as follows:
% x = unifrnd(log(1e9),log(1e11),[1e7,1]);
% [muhat,sigmahat] = normfit(x);%  
params.p.K = prior2_loggaussian('mu', 23.025, 's2', 1.33^2);
params.t.K = trans_log;
params.p.CVR = prior2_loggaussian('mu', log(0.9), 's2', 0.2);
params.t.CVR = trans_log;
params.p.R_autoc = prior2_unif('minp', minRautoc, 'maxp', maxRautoc);
params.t.R_autoc = trans_logit('minp', minRautoc, 'maxp', maxRautoc);


% Observation model, priors 
% ======================================================================= %

% variance parameter for observation model
% params.p.obs_sigma2 = prior2_sinvchi2('nu', 10, 's2', 2000^2);    
params.p.obs_sigma2 = prior2_sinvchi2('nu', 15, 's2', 1000^2);    
params.t.obs_sigma2 = trans_log;

% Length-weigth relationship parameters, weight = aw*length^bw
%   the mu values are the estimated values for log(aw) and bw given in Wieland 2002
params.p.logaw = prior2_gaussian('mu', -7.6361, 's2', 0.01);   
params.t.logaw = trans_fixed;
params.p.bw = prior2_gaussian('mu', 3.0576, 's2', 0.001); 
params.t.bw = trans_fixed;

% scaling factor (catchability) for survey index 

% survey1
params.p.q1 = prior2_loggaussian('mu', log(0.173), 's2', 0.3);
params.t.q1 = trans_log;

% survey2 = year 2003 omitted

% survey3
params.p.q3 = prior2_loggaussian('mu', log(0.173), 's2', 0.3);
params.t.q3 = trans_log;

% survey4
params.p.q4 = prior2_loggaussian('mu', log(0.173), 's2', 0.3);
params.t.q4 = trans_log;

% see Hvingel, 2013 a further discussion about this catchability:
% q is given a lognormal distribution with a median of 0.173 and a
% variance of 0.3 => sigma2 = 1.3098 (does that large a sigma2  
% make sense(?), use sigma2 = 0.3 here instead)

% Survey selectivity parameters
params.p.SL50 = prior2_loggaussian('mu',log(18), 's2', 0.005);               % length at 50% selectivity. 
params.t.SL50 = trans_log;
params.p.sigma2Ssel = prior2_sinvchi2('nu', 10, 's2', 50^2);                   % softness of the selection curve. 
params.t.sigma2Ssel = trans_log;


% Initialize the parameter vectors in the params structure
% ========================================================
fnames = fieldnames(params.p);
for i1=1:length(fnames)
    params.(fnames{i1}) = [];
end

% sample params from prior
if nargin==0
    nsamp = 1;
end
params = sample_params_from_prior(params,nsamp);

end

function state =...
   state_forward_sample(t, state, params, model, what_to_do)

I = model.I;
FL = model.FL;

% Import current state
% ====================
N = state.N;
LFmax = state.LFmax;
LR = state.LR;
nC = state.nC;
nG = state.nG;
piS = state.piS;
nR = state.nR;
Fmax = state.Fmax;

% Set fixed parameters 
% ====================
midpoints = I(1:end-1)+0.5;                     % midpoints of length classes of length 1mm           
LC = length(midpoints);                         % number of length intervals

% % Growth
minLinf = 26.5; maxLinf = 27.5;                 % The estimate given in the Pandalus stock annex is 27.1
minsdLinf = 0.01; maxsdLinf = 1;                                                                                    
minGk = 0.4; maxGk = 0.5;                       % The estimate for Gk given in the Pandalus stock annex is 0.44 

% Egg production
% NB! Pandalus is a protandric hermaphrodite, being male in the first two years (~<19mm) and female thereafter
matLength = 19;                                 % length (CL) at which maturity is reached; We assume all shrimp are mature females at 19mm and after 
mat = double(midpoints-matLength > 0);          % Maturity of a length class; 1 for mature, 0 otherwise
sex = mat;                                      % proportion of females for each length class (redundant) 


% Import parameters (with possible transformations)
% =================================================

% Growth
Linf = minLinf+(maxLinf-minLinf)*logitinv(params.logitLinf);                            
Linfsigma = minsdLinf+(maxsdLinf-minsdLinf)*logitinv(params.logitLinfsigma);            
Gk = minGk+(maxGk-minGk)*logitinv(params.logitGk);
t0 = 0;

% Natural mortality
Mconstant = params.Mconstant;

% Fishing mortality
F_autoc = params.F_autoc;                                                               % autoregression coefficient
CVF = params.CVF;
Fconstant = params.Fconstant;
L50 = params.L50;
betaGsel = params.betaGsel;

% Egg production
logaf = params.logaf;
bf = params.bf;

% Recruitment
K = params.K;
alpha = logitinv(params.logitalpha);
R_autoc = params.R_autoc;
CVR = params.CVR;



% ======================================================================= %
% Draw from prior, t=1
% ======================================================================= %
if t==1

    N = repmat(floor(K/LC),LC,1);                  % Uniform prior for distribution of size classes 
    
% Growth, t=1
% ======================================================================= %
    nG = zeros(size(N)); % initialize probability vector for population after growth

% Mortality and survival, t=1
% ======================================================================= %

% Fishing mortality, t=1
% ----------------------------------------------------------------------- %
    sigma2F = log(CVF.^2+1);
    for j1 = 1:FL
        if strcmp(what_to_do,'state_forward_sample')
            LFmax(1,1,j1) = sqrt(sigma2F(j1))*randn;
        elseif strcmp(what_to_do,'state_prior_mean')
            LFmax(1,1,j1) = muLogF;
        end
    end
    piS = zeros(size(N));
    nC = zeros(size(N));
    Fmax = zeros(size(N,2));
    
% Reproduction, t=1
% ======================================================================= %

% Egg production, t=1
% ----------------------------------------------------------------------- %
 
% Recruitment, t=1
% ----------------------------------------------------------------------- %
    sigma2R = log(CVR.^2+1);
    if strcmp(what_to_do,'state_forward_sample')
        LR = sqrt(sigma2R).*randn;                                           % stochasticity term
    elseif strcmp(what_to_do,'state_prior_mean')
        LR = 0;
    end
    nR = zeros(size(N));
    
% Update state arrays
% ===================
    state.N(:,:,:,1) = N;
    state.LFmax(:,:,:,1) = LFmax;
    state.LR(:,:,:,1) = LR;
    state.nC(:,:,:,1) = nC;
    state.nG(:,:,:,1) = nG;
    state.piS(:,:,:,1) = piS;
    state.nR(:,:,:,1) = nR;
    state.Fmax(:,:,1) = Fmax;
    return 
    
end


% ======================================================================= %
% Forward sample for t>1
% ======================================================================= %


% Growth, t>1
% ======================================================================= %
nG(:,t) = growth_expected_num2(N(:,t-1), Gk, Linfsigma, Linf, 1, midpoints, I);  % population state after growth


% Mortality and survival, t>1
% ======================================================================= %

% Natural mortality, t>1
% ----------------------------------------------------------------------- %
M = Mconstant;

% Fishing mortality, t>1
% ----------------------------------------------------------------------- %
sigma2F = log(CVF.^2+1);
mueF = Fconstant;                                                           % NB! In the manuscript, \mu_F = log(mueF)-0.5.*sigma2F
for j1 = 1:FL 
    % randomly varying F for each fleet
    if strcmp(what_to_do,'state_forward_sample')
        LFmax(1,t,j1) = F_autoc(j1)*LFmax(1,t-1,j1) +...                    % NB!LFmax is denoted \xi in the manuscript
            sqrt((1-F_autoc(j1)^2)*sigma2F(j1))*randn;  
    elseif strcmp(what_to_do,'state_prior_mean')
        LFmax(1,t,j1) = F_autoc(j1)*LFmax(1,t-1,j1);
    end
    Fmax(t,j1) = exp(log(mueF(j1))-0.5.*sigma2F(j1) + LFmax(:,t,j1));     % maximum fishing mortality for fleet j1
    Gsel = logitinv(-betaGsel(1,1,j1).*(L50(1,1,j1)-midpoints));    % gear selectivity

end

F = Fmax(t).*Gsel;

% Total mortality, t>1
% ----------------------------------------------------------------------- %
Z = M+F;

%Probabilities of mortality and survival
piS(:,t) = exp(-Z);                % probability of an individual in a given class to survive
piC = F./Z.*(1-exp(-Z));      % probability of an individual in a given class to get caught

% LLN => use expected numbers instead of draws from Bin(nG(i), piS(i)) for each class i
nS = nG(:,t).*piS(:,t);                                                               % no. of survivors per class 
nC(:,t) = nG(:,t).*piC;                                                          % no. of caught individuals per class


% Reproduction, t>1
% ======================================================================= %

% Egg production, t>1
% ----------------------------------------------------------------------- %
fec = exp(logaf+bf*log(midpoints));                                         % avg no. eggs produced by one female in each size class
eggs = sex.*mat.*fec;                                                       % avg no. eggs produced by one individual(accounting for sex and maturity) in each size class 
Eggs = N(:,t-1)'*eggs;                                                      % expected no. of eggs produced at time t, given the current population state

% Recruitment, t>1
% ----------------------------------------------------------------------- %
% \alpha: slope at the origin for the stock-recruitment function f(E)=p(E)*E, 
%   where p(E) is given by the Beverton-Holt model, \alpha \in (0,1)
% K: maximum no. of recruits

sigma2R = log(CVR.^2+1);                                                    % precision 
ER = log(Eggs) + log(K) - log(K/alpha+Eggs);                                % log of expected no. of surviving eggs (based on Beverton-Holt) relative to the initial population size 
if strcmp(what_to_do,'state_forward_sample')
    LR(1,t) = R_autoc(j1)*LR(1,t-1,j1) +...
        sqrt((1-R_autoc(j1)^2)*sigma2R(j1))*randn;
elseif strcmp(what_to_do,'state_prior_mean')
    LR(1,t) = R_autoc(j1)*LR(1,t-1,j1);
end

phiR = growth_recruits2(Gk, Linfsigma, Linf, 1, t0, midpoints, I);          % distribution of proportions of recruits over length classes
nR(:,t) = exp(ER-0.5*sigma2R + LR(1,t)).*phiR;                                   % distribution of number of recruits over length classes


% Update state arrays
% ======================================================================= %

% add survivors and recruits to obtain the state of the population
N(:,t) = nS+nR(:,t);

state.N = N;
state.LFmax = LFmax;
state.LR = LR;
state.nC = nC;
state.nG = nG;
state.piS = piS;
state.nR = nR;
state.Fmax = Fmax;
end

function lp = state_logprior(t, state, params, model, varargin)

T = model.T;
FL = model.FL;

if t>T+1
   lp=0;
   return
end

if FL >1
    error('multiple fleets not implemented yet')
end

% Import current state
% ====================
LFmax = state.LFmax;
LR = state.LR;


% Import parameters
% ====================

% Fishing mortality
F_autoc = params.F_autoc;                                                   % autoregression coefficient
CVF = params.CVF;

% Recruitment
R_autoc = params.R_autoc;
CVR = params.CVR;

% ================================================================== %
% State log prior density, t=1
% ================================================================== %
if t==1
    % LFmax
    sigma2F = log(CVF.^2+1);                                                
    lp = log(normpdf(LFmax(t), 0, sqrt(sigma2F)));

    % LR
    sigma2R = log(CVR.^2+1);                                                
    lp = lp + log(normpdf(LR(t),0,sqrt(sigma2R)));                          
    return
end


% ================================================================== %
% Forward step state log density, t>1
% ================================================================== %

%LFmax(t)|LFmax(t-1)
sigma2F = log(CVF.^2+1);
lp = log(normpdf(LFmax(t),F_autoc*LFmax(1,t-1),sqrt((1-F_autoc^2)*sigma2F)));  

%LR(t)|LR(t-1)
sigma2R = log(CVR.^2+1);                                                     
lp = lp +...
    log(normpdf(LR(t),R_autoc*LR(1,t-1),sqrt((1-R_autoc^2)*sigma2R)));
end


function ll = log_likelihood(t, state, params, model, what_to_do,... 
    totC, deltaS, surveyBiomass)

I = model.I;
FL = model.FL;
midpoints = (I(2:end)+I(1:end-1))./2;

% Import parameters
% =================
logaw = params.logaw;
bw = params.bw;
obs_sigma2 = params.obs_sigma2;
sigma2Ssel = params.sigma2Ssel;
SL50 = params.SL50;
q1 = params.q1;
q3 = params.q3;
q4 = params.q4;

% Import current state
% ====================
nC = state.nC;
nG = state.nG;
piS = state.piS;

switch what_to_do
    case 'log_likelihood'

        % If multiple fleets are to be used fix the below
        if FL >1
            error('multiple fleets not implemented yet')
        end

        
        % Observation models for catch data
        % =================================
        

       
        % total weight of catch 
        % ---------------------
        w = exp(logaw + bw*log(midpoints));                                 % mean weight of individuals per length class
        % multiply by 1000 to get the true no. of individuals,... 
        % divide by 1e6 get tons instead of g => divide by 1000
        catch_mu = nC(:,t)'*w./1000;                                        % expected weight of total estimated catch in tons
        if ~isnan(totC(t))
            wc_ll = norm_lpdf(totC(t),catch_mu,obs_sigma2);     
        else
            wc_ll = 0;
        end

        
        % Observation models for survey data
        % ==================================


        % biomass observed in survey
        % --------------------------
        % multiply by 1000 to get the true no. of individuals,... 
        % divide by 1e6 get tons instead of g => divide by 1000
        if ~isnan(surveyBiomass(t))
            Ssel = normcdf((midpoints-SL50)./sqrt(sigma2Ssel))./normcdf((31.5-SL50)./sqrt(sigma2Ssel));     % survey selectivity
            nSurv = Ssel.*nG(:,t).*(piS(:,t).^deltaS(t));                                                   % expected no. of individuals counted in the survey, up to a constant q, see observation model below
            surv_mu = nSurv'*w./1000;                                           % expected biomass of survey sample in tons
            if t>=2 && t<=20
                ws_ll = norm_lpdf(surveyBiomass(t),q1.*surv_mu,obs_sigma2); 
%           elseif t==21
%               s_ll = norm_lpdf(surveyBiomass(t),q2.*pop_mu,obs_sigma2); 
            elseif t>=22 && t<=23
                ws_ll = norm_lpdf(surveyBiomass(t),q3.*surv_mu,obs_sigma2); 
            elseif t>=24 && t<=31;
                ws_ll = norm_lpdf(surveyBiomass(t),q4.*surv_mu,obs_sigma2); 
            end
        else
            ws_ll = 0;
        end

        
        % Total log-likelihood
        % ====================
        ll = wc_ll + ws_ll;

    case 'sample'
        error('sample not implemented yet')
end

end

function state = update_state(t, state, state_new, params, model)

I = model.I;

% Import current state
% ====================
N = state.N;
LFmax = state.LFmax;
LR = state.LR;
nC = state.nC;
nG = state.nG;
piS = state.piS;
nR = state.nR;
Fmax = state.Fmax;

% Set fixed parameters 
% ====================
midpoints = I(1:end-1)+0.5;                     % midpoints of length classes of length 1mm           

% % Growth
minLinf = 26.5; maxLinf = 27.5;                 % The estimate given in the Pandalus stock annex is 27.1
minsdLinf = 0.01; maxsdLinf = 1;                                                                                    
minGk = 0.4; maxGk = 0.5;                       % The estimate for Gk given in the Pandalus stock annex is 0.44 

% Egg production
% NB! Pandalus is a protandric hermaphrodite, being male in the first two years (~<19mm) and female therafter
matLength = 19;                                 % length (CL) at which maturity is reached; We assume all shrimp are mature females at 19mm and after 
mat = double(midpoints-matLength > 0);          % Maturity of a length class; 1 for mature, 0 otherwise
sex = mat;                                      % proportion of females for each length class (redundant) 


% Import parameters (with possible transformations)
% =================================================

% Growth
Linf = minLinf+(maxLinf-minLinf)*logitinv(params.logitLinf);                            
Linfsigma = minsdLinf+(maxsdLinf-minsdLinf)*logitinv(params.logitLinfsigma);            
Gk = minGk+(maxGk-minGk)*logitinv(params.logitGk);
t0 = 0;

% Natural mortality
Mconstant = params.Mconstant;

% Fishing mortality
CVF = params.CVF;
Fconstant = params.Fconstant;
L50 = params.L50;
betaGsel = params.betaGsel;

% Egg production
logaf = params.logaf;
bf = params.bf;

% Recruitment
K = params.K;
alpha = logitinv(params.logitalpha);
CVR = params.CVR;


% ======================================================================= %
% Update state variables for t>1
% ======================================================================= %


% Growth, t>1
% ======================================================================= %
nG(:,t) = growth_expected_num2(N(:,t-1), Gk, Linfsigma, Linf, 1, midpoints, I);  % population state after growth


% Mortality and survival, t>1
% ======================================================================= %

% Natural mortality, t>1
% ----------------------------------------------------------------------- %
M = Mconstant;

% Fishing mortality, t>1
% ----------------------------------------------------------------------- %
sigma2F = log(CVF.^2+1);
mueF = Fconstant;                                                           
LFmax(t) = state_new.LFmax(t);  
Fmax = exp(log(mueF)-0.5.*sigma2F + LFmax);                              % maximum fishing mortality for fleet j1
Gsel = logitinv(-betaGsel.*(L50-midpoints));                                % gear selectivity
F = Fmax(t).*Gsel;

% Survival, t>1
% ----------------------------------------------------------------------- %
Z = M+F;
piS(:,t) = exp(-Z);                                                              % probability of an individual in a given class to survive
piC = F./Z.*(1-exp(-Z));                                                    % probability of an individual in a given class to get caught
nS = nG(:,t).*piS(:,t);                                                               % no. of survivors per class 
nC(:,t) = nG(:,t).*piC;                                                          % no. of survivors per class 

% Reproduction, t>1
% ======================================================================= %

% Egg production, t>1
% ----------------------------------------------------------------------- %
fec = exp(logaf+bf*log(midpoints));                                         % avg no. eggs produced by one female in each size class
eggs = sex.*mat.*fec;                                                       % avg no. eggs produced by one individual(accounting for sex and maturity) in each size class 
Eggs = N(:,t-1)'*eggs;                                                      % expected no. of eggs produced at time t, given the current population state

% Recruitment, t>1
% ----------------------------------------------------------------------- %
% \alpha: slope at the origin for the stock-recruitment function f(E)=p(E)*E, 
%   where p(E) is given by the Beverton-Holt model, \alpha \in (0,1)
% K: maximum no. of recruits relative to the initial population size

sigma2R = log(CVR.^2+1);                                                    % precision 
ER = log(Eggs) + log(K) - log(K/alpha+Eggs);                                % log of expected no. of surviving eggs (based on Beverton-Holt) relative to the initial population size 
LR(t) = state_new.LR(t);

phiR = growth_recruits2(Gk, Linfsigma, Linf, 1, t0, midpoints, I);          % distribution of proportions of recruits over length classes
nR(:,t) = exp(ER-0.5*sigma2R + LR(t)).*phiR;                                     % distribution of number of recruits over length classes

% Update state arrays
% ======================================================================= %

% add survivors and recruits to obtain the state of the population
N(:,t) = nS+nR(:,t);

state.N = N;
state.LFmax = LFmax;
state.LR = LR;
state.nC = nC;
state.nG = nG;
state.piS = piS;
state.nR = nR;
state.Fmax = Fmax;
end


% ================================================================
% ================================================================
% Help functions are collected below
% ================================================================
% ================================================================

% sample params from the prior
% ============================================================
function params = sample_params_from_prior(params, nsamp)
% ============================================================
for s1=1:nsamp
    % Draw samples from priors
    % ========================
    fnames = fieldnames(params.p);
    for i1=1:length(fnames)
        params.(fnames{i1})(:,s1) = params.p.(fnames{i1}).fh.rnd(params.p.(fnames{i1}));
    end
end
end

% Pack the parameters
% ===============================
function w = pack_params(params)
% ===============================
fnames = fieldnames(params.p);
for i1=1:length(fnames)
    w(i1,:) = params.t.(fnames{i1}).fh.trans(params.(fnames{i1}),params.t.(fnames{i1}));
end
end

% ================================================
function params = unpack_params(params_vec,params)
% ================================================
if isempty(params_vec)
    error('provide the params_vec')
end

fnames = fieldnames(params.p);
for i1=1:length(fnames)
    params.(fnames{i1}) = params.t.(fnames{i1}).fh.invtrans(params_vec(i1,:),params.t.(fnames{i1}));
end
end

% ===========================================
function params = take_nth_params(params,nth)
% ===========================================
fnames = fieldnames(params.p);
for i1=1:length(fnames)
    params.(fnames{i1}) = params.(fnames{i1})(:,nth);
end
end

% ========================================
function state = take_nth_state(state,nth)
% ========================================
fnames = fieldnames(state);
for i1=1:length(fnames)
    state.(fnames{i1}) = state.(fnames{i1})(:,:,:,nth);
end
end

function lp = param_logprior(params)
% ==================================
fnames = fieldnames(params.p);
lp = 0;
for i1=1:length(fnames)
    tx = params.(fnames{i1});
    p = params.p.(fnames{i1});
    t = params.t.(fnames{i1});
    lp = lp + p.fh.lp(tx,p) + t.fh.ljit(t.fh.trans(tx,t),t);
    %lp = lp + p.fh.lp(tx,p);
end
end



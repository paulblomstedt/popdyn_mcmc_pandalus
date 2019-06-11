% Copyright (c) 2014 Paul Blomstedt
% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.


%  All variables are stored in matrices for which dimensions are:
%  [LC, T, FL, S], where
%  LC = length classes
%  T  = timesteps (=years*seasons) 
%  FL = fleets
%  S = Monte Carlo samples


% Load and process the data
% =========================

% these data were provided by Mats >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
totC = load('data/tot_tons.txt');    % total yearly catch in tons
totC = [NaN*ones(1,5) totC NaN NaN]; % NB! the 1st NaN is for time step "0" and the last one is for step T+1
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%

% Estimated biomass indices from the Norwegian shrimp survey in ICES >>>>>%
% Divs. IIIa and IVa east 
survey1 = load('data/survey1.txt');
%survey2 = load('data/survey2.txt');
survey3 = load('data/survey3.txt');
survey4 = load('data/survey4.txt');
surveyBiomass = [NaN survey1 NaN survey3 survey4 NaN];
% reference:
%Søvik, G. and Thangstad, T. 2013. Results of the Norwegian Bottom Trawl 
%Survey for Northern Shrimp (Pandalus borealis) in Skagerrak and the 
%Norwegian Deep (ICES Divisions IIIa and IVa east) in 2013, 1970-2013. 
%NAFO SCR Doc. 13/071. 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%

% Fraction of the year indicating when the survey was conducted >>>>>>>>>>%
deltaS = [NaN 0.8 0.8 0.8 0.8 0.8452055 0.7958904 0.7972603 0.8041096...    % prior step, 1984:1991
    0.8041096 0.7945205 0.8287671 0.8041096 0.8027397 0.8 0.7890411...      % 1992:1998
    0.8068493 0.8 0.8027397 0.8041096 NaN 0.3767123 0.4164384 0.1041096...  % 1999:2006
    0.1191781 0.1150685 0.0890411 0.05890411 0.04794521 0.05 0.05 NaN];     % 2007:2013, step T+1
% source: stockassessment.org
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%


I = (8:32)';    % Length classes in the model in mm
T = 30+1;       % time steps in the model (30 years + 1 step for the prior)


% Initialize the model
model = model_Pandalus('T', T, 'lengthIntervals',I);


%% Run MCMC inference

% some settings

rec = [];

[~, params] = model.fh.initialize(1); 
opt.nsamp = 1000;
opt.repeat = 10;
opt.state_repeat = 10;
opt.param_repeat = 10;
opt.prop_std_s  = 0.1;
opt.prop_std = 0.01*eye(length(fieldnames(params.p)));
opt.display = 2;
opt.save = 1;

rng default % set random number generator to default settings
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 1
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.02/(sum(sum(tmp))/size(w,2)); 
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.15;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 2
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.05/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.3;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 3
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.07/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.5;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 4
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.09/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.55;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 5
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.1/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.55;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 6
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.12/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.57;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% tuning phase 7
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.14/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.6;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc


%% tuning phase 8
opt.prop_std_s  = 0.62;
opt.nsamp = 1000;
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc

%% Final run
w= rec.p';
tmp = w-ones(size(w,1),1)*mean(w);    
Sigma = tmp'*tmp./(size(w,1)-1);
tmp = diag(diag(1*chol(2.38.^2.*Sigma./size(w,2),'lower')'));
c = 0.15/(sum(sum(tmp))/size(w,2));
opt.prop_std = diag(diag(c*chol(2.38.^2.*Sigma./size(w,2),'lower')'));   % Calculate the sample variance
opt.prop_std_s  = 0.63;
opt.nsamp = 11000; % 10000 iterations + 1000 for burn-in
rng default % set random number generator to default settings
tic
[state, params, rec] = popdyn_mc2_pandalus(model, opt, rec,...
    totC, deltaS, surveyBiomass);
toc


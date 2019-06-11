function w = sample_pack(samp,w)
% pak/unpak function for state structure. 
%
% w = sample_pack(samp) return a matrix of size states x nsamples were all
% the (MCMC) sampled states are packed.
% 
% samp = sample_pack(samp,w) unpacks the packed samples from vector'
%
%  Note! unlike the packing utilities for parameters this does not
%  transform the states anyhow.

% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% pak parameters into vector
if nargin < 2
    fnames = fieldnames(samp);
    
    w=[];
    ii1=0;
    for i1=1:length(fnames)
        %ff = getfield(samp, fnames{i1});
        ff = samp.(fnames{i1});
        np = numel(ff(:,:,:,1));
        ind = ii1+1:ii1+np;
        nsamp = size(ff,4);
        w(ind,:) = reshape(ff, np,nsamp);
        ii1=ii1+np;
    end
% unpak parameters from vector
else
    nsamp = size(w,2);
    samp2 = struct;    
    fnames = fieldnames(samp);
    for i1=1:length(fnames)
        %ff = getfield(samp, fnames{i1});
        ff = samp.(fnames{i1});
        sf = size(ff);
        np = numel(ff(:,:,:,1));
        if length(sf)<3
            %samp2 = setfield(samp2, fnames{i1}, reshape(w(1:np,:),[sf(1:2) 1 nsamp]) );
            samp2.(fnames{i1}) = reshape(w(1:np,:),[sf(1:2) 1 nsamp]) ;
        else
            %samp2 = setfield(samp2, fnames{i1}, reshape(w(1:np,:),[sf(1:3) nsamp]) );
            samp2.(fnames{i1}) = reshape(w(1:np,:),[sf(1:3) nsamp]);
        end
        w(1:np,:)=[];
    end
    w=samp2;
end
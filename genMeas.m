function [mu] = genMeas(RESOL,INNER,X,epsilon,varargin)
% This code computes smoothed spectral measures of ODEs on real line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% RESOL: function handle for resolvent @(z) R_L(z)f
% INNER: function handle for inner product with f for a vector u in H (computed by RESOL) @(u) <u,f>
% X: evaluation points
% epsilon: distance to real axis from where we evaluate the resolvent

% OPTIONAL LABELLED INPUTS
% PoleType: where to place poles, default is equispaced
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% order: order of kernel used, default is 2

% OUTPUTS
% mu: smoothed measure at points X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));
addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,varargin{:})

%% Compute the higher order kernel
[pts,alpha]=rational_kernel(p.Results.Order,p.Results.PoleType);

%% Compute smoothed measure at evaluation points
mu=zeros(length(X),1);
pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));
warning('off','all')
if p.Results.Parallel=="off"
    for i=1:length(X)
        mu(i)=0;
        for ii=1:p.Results.Order
            u=alpha(ii)*RESOL(X(i)-pts(ii)*epsilon);
            mu(i)=mu(i)-imag(INNER(u))/pi;
        end
        parfor_progress(pf);
    end
else
    parfor i=1:length(X)
        mu(i)=0;
        for ii=1:p.Results.Order
            u=alpha(ii)*RESOL(X(i)-pts(ii)*epsilon);
            mu(i)=mu(i)-imag(INNER(u))/pi;
        end
        parfor_progress(pf);
    end
end
warning('on','all')
end

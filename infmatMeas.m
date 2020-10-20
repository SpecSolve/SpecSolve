function [mu] = infmatMeas(A,b,X,epsilon,varargin)
% This code computes smoothed spectral measures of finite, possibly
% rectangular matrix A. To apply this to infinite discrete matrices,
% rectangular truncations must be used so that A is f(n)\times n, where
% f describes the sparsity structure of off diagonal decay.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% A: finite rectangular matrix (sparse or full)
% b: vector in \mathbb{C}^{f(n)}
% x: vector of points where we compute smoothed measures along real line
% epsilon: the smoothing parameter

% OPTIONAL LABELLED INPUTS
% PoleType: where to place poles, default is equispaced
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% order: order of kernel used, default is 2

% OUTPUTS
% mu: smoothed measure at points x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'A',@isnumeric);
addRequired(p,'b',@isnumeric);
addRequired(p,'X',@isnumeric);
addRequired(p,'epsilon',@isnumeric);

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,A,b,X,epsilon,varargin{:})

%% first compute the poles and residues
[poles,res]=rational_kernel(p.Results.Order,p.Results.PoleType);

X=X(:);
L=length(X);
mu=zeros(size(X));
pf = parfor_progress(L);
pfcleanup = onCleanup(@() delete(pf));

if issparse(A)
    if p.Results.Parallel=="off"
        for j=1:L
            a=zeros(size(A,2),1);
            for mm=1:p.Results.Order
                a=a+res(mm)*((A-(X(j)-epsilon*poles(mm))*speye(size(A)))\b);
            end
            mu(j)=-imag(b(1:size(A,2))'*a)/pi;
            parfor_progress(pf);
        end
    else
        parfor j=1:L
            a=zeros(size(A,2),1);
            for mm=1:p.Results.Order
                a=a+res(mm)*((A-(X(j)-epsilon*poles(mm))*speye(size(A)))\b);
            end
            mu(j)=-imag(b(1:size(A,2))'*a)/pi;
            parfor_progress(pf);
        end
    end
else
    if p.Results.Parallel=="off"
        for j=1:L
            a=zeros(size(A,2),1);
            for mm=1:p.Results.Order
                a=a+res(mm)*((A-(X(j)-epsilon*poles(mm))*eye(size(A)))\b);
            end
            mu(j)=-imag(b(1:size(A,2))'*a)/pi;
            parfor_progress(pf);
        end
    else
        parfor j=1:L
            a=zeros(size(A,2),1);
            for mm=1:p.Results.Order
                a=a+res(mm)*((A-(X(j)-epsilon*poles(mm))*eye(size(A)))\b);
            end
            mu(j)=-imag(b(1:size(A,2))'*a)/pi;
            parfor_progress(pf);
        end
    end
end

end


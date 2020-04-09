function [mu] = intMeas(coeffs,f,X,epsilon,N,varargin)
% This code computes smoothed spectral measures of an integral operator 
% with form Lu(x)=a(x)u(x)+int_{-1}^1 G(x,y)u(y)dy,  x \in [-1,1].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% coeffs: cell of op. coefficients given as function handles, i.e. 
%         coeffs={@(x) a(x), @(x,y) G(x,y)}
% f: function handle for function we compute measure wrt, i.e. mu_f
% X: evaluation points
% order: order of kernel used
% epsilon: distance to real axis from where we evaluate the resolvent
% N: size of N x N discretization from Chebyshev collocation scheme

% OPTIONAL LABELLED INPUTS
% PoleType: where to place poles, default is equispaced
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% order: order of kernel used, default is 2

% OUTPUTS
% mu: smoothed measure wrt f at points X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'coeffs',@iscell);
addRequired(p,'f',@(x) isa(x,'function_handle'));
addRequired(p,'X',@isnumeric);
addRequired(p,'epsilon',@isnumeric);
addRequired(p,'N',@isnumeric);

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,coeffs,f,X,epsilon,N,varargin{:})

%% Multiplication by a(x) and kernel G(x,y) with rank k representation
a=coeffs{1};
G=chebfun2(coeffs{2}); [Gc,Gd,Gr]=cdr(G);

%% Construct discretizations
[ccn,ccw]=chebpts(N,1);

%Discretized integral operator with low-rank kernel and right-hand side
M_a=a(ccn);           %multiplicative term
M_Gr=ccw.*Gr(ccn)';   %kernel row basis (y) with quad weights
M_Gc=Gc(ccn);         %kernel column basis (x)
f_vals=f(ccn);        %Right hand side
Gd_inv=diag(1./diag(Gd));   %invert Gd for Woodbury solve

%% Compute smoothed measure at evaluation points
[poles,res]=rational_kernel(p.Results.Order,p.Results.PoleType);

%Solve diagonal plus low-rank with Woodbury formula
mu=zeros(size(X));
pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));
warning('off','all')
if p.Results.Parallel=="off"
    for n=1:length(X)
        u=zeros(N,1);
        for j=1:p.Results.Order
            z=X(n)+epsilon*poles(j);
            u1=f_vals./(M_a-z);
            u2=M_Gc./(M_a-z);
            u=u+res(j)*(u1-u2*((Gd_inv+M_Gr*u2)\(M_Gr*u1)));
        end

        %Compute inner product with Clenshaw-Curtis quadrature
        mu(n)=imag(ccw*(u.*conj(f_vals)))/pi;
        parfor_progress(pf);
    end
else
    parfor n=1:length(X)
        u=zeros(N,1);
        for j=1:p.Results.Order
            z=X(n)+epsilon*poles(j);
            u1=f_vals./(M_a-z);
            u2=M_Gc./(M_a-z);
            u=u+res(j)*(u1-u2*((Gd_inv+M_Gr*u2)\(M_Gr*u1)));
        end

        %Compute inner product with Clenshaw-Curtis quadrature
        mu(n)=imag(ccw*(u.*conj(f_vals)))/pi;
        parfor_progress(pf);
    end
end
warning('on','all')
end
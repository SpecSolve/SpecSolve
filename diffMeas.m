function [mu] = diffMeas(a,f,X,epsilon,N,varargin)
% This code computes smoothed spectral measures of ODEs on real line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% a: cell of coefficients given as function handles so that L=a_0
% +a_1*d/dx+...+a_p(d/dx)^p
% f: function that we compute measure with respect to, i.e. mu_f
% X: evaluation points
% epsilon: distance to real axis from where we evaluate the resolvent
% N: number of Fourier modes (actually 2*N+1)

% OPTIONAL LABELLED INPUTS
% Band_max: max bandwidth for coefficients, default is 200 
% scale: scale for mapping to [-pi,pi], default is 10
% PoleType: where to place poles, default is equispaced
% Parallel: parfor (on) or normal for (off) loop, default is "off"
% order: order of kernel used, default is 2

% OUTPUTS
% mu: smoothed measure at points X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'a',@iscell);
addRequired(p,'f',@(x) isa(x,'function_handle'));
addRequired(p,'X',@isnumeric);
addRequired(p,'epsilon',@isnumeric);
addRequired(p,'N',@isnumeric);

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'Band_max',min(200,N),@(x) x>0)
addParameter(p,'Scale',10,@(x) x>0)
addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,a,f,X,epsilon,N,varargin{:})

%% Compute the higher order kernel
[pts,alpha]=rational_kernel(p.Results.Order,p.Results.PoleType);

%% Build discretized linear system for solution of shifted ODEs
% Map coefficient and right-hand-side functions to unit circle
map=@(t) p.Results.Scale*1i*(1-exp(1i*t))./(1+exp(1i*t));
singular_weight=chebfun(@(t) 0.5*exp(-1i*t).*(1+exp(1i*t)).^2,[-pi,pi],'trig')/p.Results.Scale;
F=chebfun(@(t) f(map(t)),[-pi,pi],'trig');
G=chebfun(@(t) p.Results.Scale*2*exp(1i*t).*f(map(t))./(1+exp(1i*t)).^2,[-pi,pi],'trig');

% Construct differentiation matrix
D=1i*spdiags((-N:N)',0,2*N+1,2*N+1);
weight_coeffs=trigcoeffs(singular_weight,4*N+1);
SW=sptoeplitz(weight_coeffs(2*N+1:end),weight_coeffs(2*N+1:-1:1));
D=SW*D;

% Construct variable coeff diff op matrix
for k=1:length(a)
    var_coeffs=trigcoeffs(chebfun(@(t) a{k}(map(t)),[-pi,pi],'trig'),2*p.Results.Band_max+1);
    if 2*N>p.Results.Band_max
        var_coeffs=[sparse(2*N-p.Results.Band_max,1);
                    var_coeffs;
                    sparse(2*N-p.Results.Band_max,1)];
    end
    C=sptoeplitz(var_coeffs(2*N+1:end),var_coeffs(2*N+1:-1:1));
    if k==1
        L=C; % extra "if" needed, otherwise get memory error
    else
        L=L+C*(D^(k-1));
    end
end

% Discretize right-hand-side
F_coeffs=sqrt(2*pi)*trigcoeffs(F,2*N+1);
G_coeffs=sqrt(2*pi)*trigcoeffs(G,2*N+1);

%% Compute smoothed measure at evaluation points
mu=zeros(length(X),1); I=speye(2*N+1);
pf = parfor_progress(length(X));
pfcleanup = onCleanup(@() delete(pf));
warning('off','all')
if p.Results.Parallel=="off"
    for i=1:length(X)
        u_coeffs=zeros(size(F_coeffs));
        for ii=1:p.Results.Order
            u_coeffs=u_coeffs+alpha(ii)*((L-(X(i)-pts(ii)*epsilon)*I)\F_coeffs);
        end
        mu(i)=-imag(G_coeffs'*u_coeffs)/pi;
        parfor_progress(pf);
    end
else
    parfor i=1:length(X)
        u_coeffs=zeros(size(F_coeffs));
        for ii=1:p.Results.Order
            u_coeffs=u_coeffs+alpha(ii)*((L-(X(i)-pts(ii)*epsilon)*I)\F_coeffs);
        end
        mu(i)=-imag(G_coeffs'*u_coeffs)/pi;
        parfor_progress(pf);
    end
end
warning('on','all')
end


function [mu] = diffMeas(a,f,X,epsilon,varargin)
% This code computes smoothed spectral measures of ODEs on real line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% a: cell of coefficients given as function handles so that L=a_0
% +a_1*d/dx+...+a_p(d/dx)^p
% f: function that we compute measure with respect to, i.e. mu_f
% X: evaluation points
% epsilon: distance to real axis from where we evaluate the resolvent

% OPTIONAL LABELLED INPUTS
% N: number of Fourier modes (discretization size actually 2*N+1)
% DiscMin: maximum discretization size
% DiscMin: minimum discretization size
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

validPole = {'cheb','roots','extrap','equi'};
checkPole = @(x) any(validatestring(x,validPole));
validPar = {'on','off'};
checkPar = @(x) any(validatestring(x,validPar));

addParameter(p,'N',[],@(x) x==floor(x))
addParameter(p,'DiscMin',2^10,@(x) x==floor(x))
addParameter(p,'DiscMax',2^15,@(x) x==floor(x))
addParameter(p,'Band_max',200,@(x) x>0)
addParameter(p,'Scale',10,@(x) x>0)
addParameter(p,'PoleType','equi',checkPole)
addParameter(p,'Parallel','off',checkPar)
addParameter(p,'Order',2,@(x) x==floor(x))

p.CaseSensitive = false;
parse(p,a,f,X,epsilon,varargin{:})

% Get eval pts and weights for rational convolution kernels
[pts,alpha]=rational_kernel(p.Results.Order,p.Results.PoleType);

% Set up chebfuns for variable coeffs and rhs
map=@(t) p.Results.Scale*1i*(1-exp(1i*t))./(1+exp(1i*t));
singular_weight=chebfun(@(t) 0.5*exp(-1i*t).*(1+exp(1i*t)).^2,[-pi,pi],'trig')/p.Results.Scale;
F=chebfun(@(t) f(map(t)),[-pi,pi],'trig');
G=chebfun(@(t) p.Results.Scale*2*exp(1i*t).*f(map(t))./(1+exp(1i*t)).^2,[-pi,pi],'trig');
var_coeffs=zeros(2*p.Results.Band_max+1,length(a));
for k=1:length(a)
    var_coeffs(:,k)=trigcoeffs(chebfun(@(t) a{k}(map(t)),[-pi,pi],'trig'),2*p.Results.Band_max+1);
end

if isempty(p.Results.N)
    % Adaptive resolvent evaluations on grid
    N=p.Results.DiscMin;                           %init discSize
    mu=zeros(size(X)); errIndx=1:length(X);
    tol=epsilon^(p.Results.Order+1);
    while(sum(errIndx)>0 && N<=p.Results.DiscMax)   %refine measure     

        % Construct differentiation matrix
        D=1i*spdiags((-N:N)',0,2*N+1,2*N+1);
        weight_coeffs=trigcoeffs(singular_weight,4*N+1);
        SW=sptoeplitz(weight_coeffs(2*N+1:end),weight_coeffs(2*N+1:-1:1));
        D=SW*D;

        % Construct variable coeff diff op matrix
        for k=1:length(a)
            if 2*N>p.Results.Band_max
                var_coeffsT=[sparse(2*N-p.Results.Band_max,1);
                            var_coeffs(:,k);
                            sparse(2*N-p.Results.Band_max,1)];
                C=sptoeplitz(var_coeffsT(2*N+1:end),var_coeffsT(2*N+1:-1:1));
            else
                C=sptoeplitz(var_coeffs(2*N+1:end,k),var_coeffs(2*N+1:-1:1,k));
            end
            if k==1
                L=C; % extra "if" needed, otherwise get memory error
            else
                L=L+C*(D^(k-1));
            end
        end

        % Shifts
        I=speye(2*N+1);

        % Discretize right-hand-side
        F_coeffs=sqrt(2*pi)*trigcoeffs(F,2*N+1);
        G_coeffs=sqrt(2*pi)*trigcoeffs(G,2*N+1);

        % Compute smoothed measure at evaluation points
        if p.Results.Parallel=="off"
            temp=mu;
            for i=1:length(X(errIndx))
                u_coeffs=zeros(size(F_coeffs));
                for ii=1:p.Results.Order
                    u_coeffs=u_coeffs+alpha(ii)*((L-(X(errIndx(i))-pts(ii)*epsilon)*I)\F_coeffs);
                end
                mu(errIndx(i))=-imag(G_coeffs'*u_coeffs)/pi;
            end
            errIndx=find(abs(mu-temp)>tol*abs(mu));
            N=2*N;
        else
            pf = parfor_progress(length(X));
            pfcleanup = onCleanup(@() delete(pf));
            warning('off','all')
            temp=mu(errIndx); muTemp=temp; XTemp=X(errIndx);
            parfor i=1:length(XTemp)
                u_coeffs=zeros(size(F_coeffs));
                for ii=1:p.Results.Order
                    u_coeffs=u_coeffs+alpha(ii)*((L-(XTemp(i)-pts(ii)*epsilon)*I)\F_coeffs);
                end
                muTemp(i)=-imag(G_coeffs'*u_coeffs)/pi;
                parfor_progress(pf);
            end
            errIndx=find(abs(muTemp-temp)>tol*abs(muTemp));
            mu(errIndx)=muTemp;
            N=2*N;
            warning('on','all')
        end
    end
else
    % Fixed N resolvent evaluations on grid
    N=p.Results.N;
    mu=zeros(size(X));

    % Construct differentiation matrix
    D=1i*spdiags((-N:N)',0,2*N+1,2*N+1);
    weight_coeffs=trigcoeffs(singular_weight,4*N+1);
    SW=sptoeplitz(weight_coeffs(2*N+1:end),weight_coeffs(2*N+1:-1:1));
    D=SW*D;

    % Construct variable coeff diff op matrix
    for k=1:length(a)
        if 2*N>p.Results.Band_max
            var_coeffsT=[sparse(2*N-p.Results.Band_max,1);
                        var_coeffs(:,k);
                        sparse(2*N-p.Results.Band_max,1)];
            C=sptoeplitz(var_coeffsT(2*N+1:end),var_coeffsT(2*N+1:-1:1));
        else
            C=sptoeplitz(var_coeffs(2*N+1:end,k),var_coeffs(2*N+1:-1:1,k));
        end
        if k==1
            L=C; % extra "if" needed, otherwise get memory error
        else
            L=L+C*(D^(k-1));
        end
    end

    % Shifts
    I=speye(2*N+1);

    % Discretize right-hand-side
    F_coeffs=sqrt(2*pi)*trigcoeffs(F,2*N+1);
    G_coeffs=sqrt(2*pi)*trigcoeffs(G,2*N+1);

    % Compute smoothed measure at evaluation points
    if p.Results.Parallel=="off"
        pf = parfor_progress(length(X));
        pfcleanup = onCleanup(@() delete(pf));
        warning('off','all')
        for i=1:length(X)
            u_coeffs=zeros(size(F_coeffs));
            for ii=1:p.Results.Order
                u_coeffs=u_coeffs+alpha(ii)*((L-(X(i)-pts(ii)*epsilon)*I)\F_coeffs);
            end
            mu(i)=-imag(G_coeffs'*u_coeffs)/pi;
            parfor_progress(pf);
        end
    else
        pf = parfor_progress(length(X));
        pfcleanup = onCleanup(@() delete(pf));
        warning('off','all')
        parfor i=1:length(X)
            u_coeffs=zeros(size(F_coeffs));
            for ii=1:p.Results.Order
                u_coeffs=u_coeffs+alpha(ii)*((L-(X(i)-pts(ii)*epsilon)*I)\F_coeffs);
            end
            mu(i)=-imag(G_coeffs'*u_coeffs)/pi;
            parfor_progress(pf);
        end
        warning('on','all')
    end
end
end
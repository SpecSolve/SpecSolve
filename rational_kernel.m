function [poles,res] = rational_kernel(m,type)
%%rational_kernel 
% Inputs: 
%         - order of kernel 'm'
%         - Pole location 'roots','cheb'
% Outputs: - poles of rational function 'poles' 
%         - residues of rational function 'res'

%Poles
if type=="cheb"
    z=1i+chebpts(m,1);  %Chebyshev points
    %z=1i+(2*(1:m)/(m+1)-1);%.^2.*sign((2*(1:m)/(m+1)-1));
elseif type=="roots"
    z=exp(1i*pi*((1:m)'-1/2)/m); %Roots of unity
elseif type=="extrap"
    z=1i*(0.5.^(0:m-1));
elseif type=="equi"
    z=1i+(2*(1:m)/(m+1)-1);
end

if type=="equi" && m<7 %Hard-coded kernels from table of the paper
    if m==1
        res=1;
    elseif m==2
        res=[(1+3i)/2;(1-3i)/2];
    elseif m==3
        res=[-2+1i;5;-2-1i];
    elseif m==4
        res=[(-39-65i)/24;(17+85i)/8;(17-85i)/8;(-39+65i)/24];
    elseif m==5
        res=[(15-10i)/4;(-39+13i)/2;(65)/2;(-39-13i)/2;(15+10i)/4];
    else
        res=[(725+1015i)/(192);(-2775-6475i)/(192);(1073+7511i)/(96);(1073-7511i)/(96);(-2775+6475i)/(192);(725-1015i)/(192)];
    end
else %Vandermonde matrix for poles
    V=zeros(m);
    for i=1:m 
        V(:,i)=(z).^(i-1);
    end

    %Get residues directly
    rhs=eye(m,1);
    res=transpose(V)\rhs;
end

poles=z;

end


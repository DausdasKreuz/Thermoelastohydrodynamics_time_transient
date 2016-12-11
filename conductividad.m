% Cálculo de la conductividad térmica para cada punto del fluido en el
% instante dado mediante un modelo combinado de temperatura y de presión.


function [k]=conductividad(kini,T,p,C_k1,C_k2,npf)

k=zeros(npf,1);

for f=1:npf
    k(f)=kini*(1+C_k1*p/(1+C_k2*p))*(1-1.667e-4*T(f));
end
% Cálculo de la densidad para cada punto del fluido en el instante dado
% mediante el modelo de presión y temperatura de Dowson-Higginson.


function [ro]=densidad(roini,p,T,Tini,epsilon,C_ro1,C_ro2,npf)

ro=zeros(npf,1);

for f=1:npf
    ro(f)=roini*(1+(C_ro1*p)/(1+C_ro2*p))*(1-epsilon*(T(f)-Tini));
end
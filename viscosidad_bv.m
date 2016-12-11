% Cálculo de la viscosidad a baja velocidad para cada punto del fluido en
% el instante dado mediante el modelo de Roelands.


function [mu]=viscosidad_bv(muini,p,T,Tini,alfa_mu,beta_mu,npf)

mu=zeros(npf,1);

z0=alfa_mu/(5.1e-9*(log(muini)+9.67));
s0=beta_mu*(Tini-138)/(log(muini)+9.67);
for f=1:npf
    mu(f)=muini*exp((log(muini)+9.67)*(-1+(1+5.1e-9*p)^z0*...
        ((T(f)-138)/(Tini-138))^(-s0)));
end
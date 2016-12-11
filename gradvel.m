% Cálculo del gradiente de velocidad para cada punto del fluido en el
% instante dado mediante el modelo de Ree-Eyring.


function [gdv]=gradvel(tau,tau0,eta,npf)

gdv=zeros(npf,1);

for f=1:npf
    gdv(f)=tau0/eta(f)*sinh(tau(f)/tau0);
end
% Cálculo del esfuerzo cortante para cada punto del fluido en el instante
% dado mediante el modelo de Ree-Eyring.


function [tau]=cortante(tau0,eta,gdv,npf)

tau=zeros(npf,1);
tauacumulado=0;

for f=1:npf
    tau(f)=tau0*asinh(eta(f)*gdv(f)/tau0);
    tauacumulado=tauacumulado+tau(f);
end

tau=ones(npf,1)*(tauacumulado/npf);
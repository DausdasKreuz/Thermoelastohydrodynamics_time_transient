% C�lculo de la tensi�n l�mite del lubricante durante el transitorio de
% presi�n a partir de un modelo generalizado.

function [taulim]=tau_limite(tau0,dseda,p,npt)

taulim=zeros(npt,1);

for i=1:npt
    taulim(i)=tau0+dseda*p(i);
end
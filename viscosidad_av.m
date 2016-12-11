% Cálculo de la viscosidad generalizada (a alta velocidad) para cada punto
% del fluido en el instante dado medianteel modelo de Carreau.


function [eta]=viscosidad_av(mu,gdv,G,npf,cCarr,expley)

eta=zeros(npf,1);

for f=1:npf
    eta(f)=mu(f)*(1+(mu(f)*gdv(f)/G)^cCarr)^((expley-1)/cCarr);
end
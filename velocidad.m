% Cálculo de la velocidad de los puntos del fluido en cada momento.

function [v]=velocidad(gdv,hz,npf)

v=zeros(npf,1);

for f=2:npf
    v(f)=gdv(f)*hz+v(f-1);
end
% Cálculo del calor específico para cada punto del fluido en el instante
% dado mediante el modelo de Larsson-Andersson


function [cp]=calor_esp(cpini,roini,Tini,ro,T,p,k_cp1,k_cp2,beta_cpini,...
    b_cp1,b_cp2,npf)

cp=zeros(npf,1);

beta_cp=beta_cpini*(1+b_cp1*(p*1e-9)+b_cp2*(p*1e-9)^2);
for f=1:npf
    cp(f)=(roini*cpini*(1+(k_cp1*p/(1+k_cp2*p)))*(1+beta_cp*...
        (T(f)-Tini)))/ro(f);
end
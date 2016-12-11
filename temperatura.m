% Cálculo de la temperatura en todo el fluido considerado para el instante
% dado.  Para la obtención de esta magnitud se crea un sistema de
% ecuaciones en el que intervienen las temperaturas de cada uno de los
% puntos del fluido.


function [T]=temperatura(Tm1,Tinf,Tsup,ro,rom1,cp,ht,k,hz,tau,gdv,...
    epsilon,dp,npf)

dro=(ro-rom1)/ht;   %Cálculo de la derivada de la densidad.

% Cálculo de los coeficientes de la matriz del sistema de ecuaciones.
S=zeros(npf,npf);       %matriz del sistema de ecuaciones
cT=zeros(npf,1);        %coef del término T
cdT=zeros(npf,1);       %coef de los términos T+-j
cindep=zeros(npf,1);    %términos independientes de las ec's del sistema

%Vectores de coeficientes de la matriz del sistema
% (derivada espacial central de segundo orden)
for f=2:npf-1
    cT(f)=(cp(f)*dro(f))+(cp(f)*ro(f)/ht)+(2*k(f)/hz^2)-(epsilon*dp);
    cdT(f)=-k(f)/hz^2;
end

% Creación de la matriz del sistema
S(1,1)=1;       %Condición de evacuación del calor (condición de Neumann)
S(npf,npf)=1;   %Condición de evacuación del calor (condición de Neumann)

for s=2:npf-1
    S(s,s-1)=cdT(f);
    S(s,s)=cT(f);
    S(s,s+1)=cdT(f);
end

% Creación del vector de términos independientes
cindep(1)=Tinf;     %(condición de Neumann)
cindep(npf)=Tsup;   %(condición de Neumann)

for f=2:npf-1       %Resto de términos independientes
    cindep(f)=(cp(f)*ro(f)*Tm1(f))/ht+tau(f)*gdv(f);
end

% Cálculo de las temperaturas: resolución del sistema --> S*T=cindep
T=S\cindep;
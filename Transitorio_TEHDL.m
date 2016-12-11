% Aplicación que obtiene los valores de las variables que intervienen en el
% comportamiento de los contactos lubricados con movimiento relativo entre
% los sólidos que lo delimitan.
%
% En esta parte se introducen los siguientes datos datos:
%   -las variables matemáticas necesarias para la aplicación
%   -el transitorio de carga a que es sometido el sistema
%   -el lubricante que funciona como medio fluido
%   -el sistema físico que simula el contacto lubricado
% 
% Finalmente se genera un vector con todos los datos y se llama a la
% función que va a efectuar los cálculos devolviendo los resultados
% finales.



% Variables matemáticas
nit=1e3;            % número de iteraciones máximo
error=1e-7;         % error límite para validaro los resultados
ndt=100;            % número de divisiones del intervalo temporal
ndf=100;            % número de divisiones en el espesor del fluido


% Datos del transitorio
tini=0;             % tiempo del instante inicial [s]
tend=0.0002;        % tiempo del instante final [s]


% Datos del lubricante
% Tipo de lubricante: PAO
lubricante=10;   % Elegir PAO: 4->PAO4; 10->PAO10; 40->PAO40; 100->PAO100
[cpini,kini,roini,muini,G,dseta,epsilonini,c_epsilon,k_cp1,k_cp2,...
    beta_cpini,b_cp1,b_cp2,C_k1,C_k2,alfa_mu,beta_mu,C_ro1,C_ro2,tau0,...
    cCarr,expley]=familia_PAO(lubricante);


% Datos del sistema
h=1e-6;             % espesor==distancia entre las placas [m]
gdvini=1e6;         % gradiente de velocidad [m/(s*m)]
pmax=2e9;           % presión máxima a que es sometido el sistema [Pa]
Tini=298;           % temperatura inicial del sistema [K]
Tinf=298;           % temperatura de evacuación inferior [K]
Tsup=298;           % temperatura de evacuación superior [K]


[T,tau,gdv,v,vrel,mur,mu,eta,ro,k,cp,epsilon,taulim,p,dp,N]=...
    Calculos_TEHDL(nit,error,ndt,ndf,tini,tend,cpini,kini,roini,muini,G,...
    dseta,epsilonini,c_epsilon,k_cp1,k_cp2,beta_cpini,b_cp1,b_cp2,C_k1,...
    C_k2,alfa_mu,beta_mu,C_ro1,C_ro2,tau0,cCarr,expley,h,gdvini,pmax,...
    Tini,Tinf,Tsup);
% Aplicaci�n que obtiene los valores de las variables que intervienen en el
% comportamiento de los contactos lubricados con movimiento relativo entre
% los s�lidos que lo delimitan.
%
% En esta parte se introducen los siguientes datos datos:
%   -las variables matem�ticas necesarias para la aplicaci�n
%   -el transitorio de carga a que es sometido el sistema
%   -el lubricante que funciona como medio fluido
%   -el sistema f�sico que simula el contacto lubricado
% 
% Finalmente se genera un vector con todos los datos y se llama a la
% funci�n que va a efectuar los c�lculos devolviendo los resultados
% finales.



% Variables matem�ticas
nit=1e3;            % n�mero de iteraciones m�ximo
error=1e-7;         % error l�mite para validaro los resultados
ndt=100;            % n�mero de divisiones del intervalo temporal
ndf=100;            % n�mero de divisiones en el espesor del fluido


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
pmax=2e9;           % presi�n m�xima a que es sometido el sistema [Pa]
Tini=298;           % temperatura inicial del sistema [K]
Tinf=298;           % temperatura de evacuaci�n inferior [K]
Tsup=298;           % temperatura de evacuaci�n superior [K]


[T,tau,gdv,v,vrel,mur,mu,eta,ro,k,cp,epsilon,taulim,p,dp,N]=...
    Calculos_TEHDL(nit,error,ndt,ndf,tini,tend,cpini,kini,roini,muini,G,...
    dseta,epsilonini,c_epsilon,k_cp1,k_cp2,beta_cpini,b_cp1,b_cp2,C_k1,...
    C_k2,alfa_mu,beta_mu,C_ro1,C_ro2,tau0,cCarr,expley,h,gdvini,pmax,...
    Tini,Tinf,Tsup);
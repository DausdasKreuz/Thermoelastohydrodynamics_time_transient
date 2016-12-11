% Fichero de datos de lubricantes pertenecientes a la familia de las
% polialfaolefinas (PAO).
% En la llamada a este archivo se ha seleccionado la PAO concreta para el
% estudio.


function [cp,k,ro,mu,G,dseta,epsilon,c_epsilon,k_cp1,k_cp2,beta_cp,...
    b_cp1,b_cp2,C_k1,C_k2,alfa_mu,beta_mu,C_ro1,C_ro2,tau0,cCarr,expley]...
    =familia_PAO(lubricante)

% Selección del lubricante concreto
if lubricante==4
    datos_PAO4;
elseif lubricante==10
    datos_PAO10;
elseif lubricante==40
    datos_PAO40;
elseif lubricante==100
    datos_PAO100;
else
    disp('No se ha seleccionado convenientemente el lubricante.')
end

% Datos generales de la familia de (PAO)
G=1.59e9;
dseta=0.0395;
epsilon=6.8e-4;
c_epsilon=1.5e-9;
k_cp1=0.41e-9;
k_cp2=1.05e-9;
beta_cp=6.5e-4;
b_cp1=2.7;
b_cp2=-1.5;
C_k1=1.40e-9;
C_k2=0.34e-9;
beta_mu=0.03;
C_ro1=0.69e-9;
C_ro2=2.55e-9;
tau0=4e6;
cCarr=2;
expley=0.2;
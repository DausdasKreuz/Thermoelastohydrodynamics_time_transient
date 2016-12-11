% Se definen los valores de presión a lo largo del tiempo en función de los
% parámetros del transitorio.

% La presión es la condición que se le impone al sistema, por lo que la
% función de la presión se define de forma arbitraria para obtener un pulso
% con cambios suaves que no presenten puntos en los que la variación de
% magnitudes sea demasiado grande.

% Para el cálculo de la derivada se emplea la expresión analítica de la
% derivada de la presión, ya que de esta forma será más precisa
% independientemente del paso que se tome.


function [p,dp]=presion(pmax,npt,ndt,ht)

p=zeros(npt,1);    % vector de la presión
dp=zeros(npt,1);    % vector derivada de la presión

% Tiempo de relajamiento antes y después = 1/4 tiempo transitorio
% Cálculo de la presión
for i=ndt*0.25:ndt*1.25
    p(i)=pmax*(sin(((i-ndt)/ndt*2*pi)+pi)+1)/2;
end

% Cálculo de la derivada de la presión
for i=ndt*0.25:ndt*1.25
    dp(i)=pmax*pi*(cos(((i-ndt)/ndt*2*pi)+pi))/2;
end
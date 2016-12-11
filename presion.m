% Se definen los valores de presi�n a lo largo del tiempo en funci�n de los
% par�metros del transitorio.

% La presi�n es la condici�n que se le impone al sistema, por lo que la
% funci�n de la presi�n se define de forma arbitraria para obtener un pulso
% con cambios suaves que no presenten puntos en los que la variaci�n de
% magnitudes sea demasiado grande.

% Para el c�lculo de la derivada se emplea la expresi�n anal�tica de la
% derivada de la presi�n, ya que de esta forma ser� m�s precisa
% independientemente del paso que se tome.


function [p,dp]=presion(pmax,npt,ndt,ht)

p=zeros(npt,1);    % vector de la presi�n
dp=zeros(npt,1);    % vector derivada de la presi�n

% Tiempo de relajamiento antes y despu�s = 1/4 tiempo transitorio
% C�lculo de la presi�n
for i=ndt*0.25:ndt*1.25
    p(i)=pmax*(sin(((i-ndt)/ndt*2*pi)+pi)+1)/2;
end

% C�lculo de la derivada de la presi�n
for i=ndt*0.25:ndt*1.25
    dp(i)=pmax*pi*(cos(((i-ndt)/ndt*2*pi)+pi))/2;
end
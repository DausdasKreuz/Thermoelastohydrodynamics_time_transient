% C�lculo del coeficiente de expansi�n t�rmica del lubricante durante el
% transitorio de presi�n a partir de un modelo generalizado.

function [epsilon]=coef_expanion(epsilonini,c_epsilon,p,npt)

epsilon=zeros(npt,1);

for i=1:npt
    epsilon(i)=epsilonini*exp(-c_epsilon*p(i));
end
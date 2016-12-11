% Cálculo del coeficiente de expansión térmica del lubricante durante el
% transitorio de presión a partir de un modelo generalizado.

function [epsilon]=coef_expanion(epsilonini,c_epsilon,p,npt)

epsilon=zeros(npt,1);

for i=1:npt
    epsilon(i)=epsilonini*exp(-c_epsilon*p(i));
end
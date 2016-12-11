% Función que calcula las principales variables que intervienen en un
% contacto lubricado. Se obtienen los valores de viscosidad (mu),
% viscosidad generalizada (eta), densidad (ro), temperatura (T), esfuerzo
% cortante (tau), gradiente de velocidad (gdv), velocidad del fluido (v),
% ceficiente de fricción (mur), coeficiente de expansión (epsilon) y
% esfuerzo cortante límite (taulim).

% El proceso se desarrolla en el tiempo, estando sujeto a una carga de
% valor variable en el tiempo, lo que dando lugar a un transitorio
% temporal. En t<t0 (i<=2) el sistema  se encuentra a temperatura "Tini".
% Desde el instante t=t0 (i=3) hasta el instante i=1/4*ndt el sistema se
% encuentra en reposo. Al llegar el instante i=(1/4+1)*ndt el transitorio
% termina quedando en reposo nuevamente el sistema.

% Durante todo el proceso se evacua calor o se calienta el lubricante a
% través de las paredes de los elementos sólidos que actuan de frontera.

% Se ha utilizado un método explicito (diferencia 'backwards' de la
% derivada temporal).

% Los cálculos se realizan tanto en el espacio (h --> "f") como en
% el tiempo (t --> "i") dando lugar a resultados en forma de superficie.



function [T,tau,gdv,v,vrel,mur,mu,eta,ro,k,cp,epsilon,taulim,p,dp,N]=...
    Calculos_TEHDL(nit,error,ndt,ndf,tini,tend,cpini,kini,roini,muini,G,...
    dseta,epsilonini,c_epsilon,k_cp1,k_cp2,beta_cpini,b_cp1,b_cp2,C_k1,...
    C_k2,alfa_mu,beta_mu,C_ro1,C_ro2,tau0,cCarr,expley,h,gdvini,pmax,...
    Tini,Tinf,Tsup)


% 1. -------------Cálculo de variables auxiliares--------------%
npt=1.5*ndt+1;     % número de puntos del intervalo
npf=ndf+1;         % número de puntos del fluido

ht=(tend-tini)/ndt;     % paso temporal [s]
hz=h/ndf;               % paso en el fluido [m]



% 2. -------------Reserva de espacio de las variables--------------%
% espacio-h, tiempo-t
cp=zeros(npf,npt);          % Matriz de calor específico (h,t)
k=zeros(npf,npt);           % Matriz de conductividad térmica (h,t)
mu=zeros(npf,npt);          % Matriz de viscosidad (h,t)
eta=zeros(npf,npt);         % Matriz de viscosidad generalizada (h,t)
ro=zeros(npf,npt);          % Matriz de densidad (h,t)
tau=zeros(npf,npt);         % Matriz de esfuerzo cortante (h,t)
gdv=zeros(npf,npt);         % Matriz de gradiente de velocidad (h,t)
v=zeros(npf,npt);           % Matriz de velocidad del fluido (h,t)
vrel=zeros(npf,npt);        % Matriz de velocidad relativa (h,t)
mur=zeros(npf,npt);         % Matriz de coeficiente de rozamiento (h,t)
T=zeros(npf,npt);           % Matriz de temperatura (h,t)
difT=(error+1)*ones(npf,1); % diferencia entre valores de iteraciones

N=zeros(npt,1);             %nº de iteraciones en cada instante



% 3. -------------Cálculo de variables físicas--------------%
% 3.1 * Variables del sistema
[p,dp]=presion(pmax,npt,ndt,ht);  % valores de presión

[taulim]=tau_limite(tau0,dseta,p,npt);  %valores de esfuerzo corteante lím.
[epsilon]=coef_expanion(epsilonini,c_epsilon,p,npt); % coef. expansión


% 3.2 * Instantes iniciales: desde i=1 hasta i=2
cp(:,1)=cpini*ones(npf,1);
k(:,1)=kini*ones(npf,1);
mu(:,1)=muini*ones(npf,1);
gdv(:,1)=gdvini*ones(npf,1);
[eta(:,1)]=viscosidad_av(mu(:,1),gdv(:,1),G,npf,cCarr,expley);
ro(:,1)=roini*ones(npf,1);
[tau(:,1)]=cortante(tau0,eta(:,1),gdv(:,1),npf);
[v(:,1)]=velocidad(gdv(:,1),hz,npf);
vrel(:,1)=v(:,1)./v(:,1);
mur(:,1)=tau(:,1)/p(1);
T(:,1)=Tini*ones(npf,1);

cp (:,2)= cp (:,1);
k  (:,2)= k  (:,1);
mu (:,2)= mu (:,1);
eta(:,2)= eta(:,1);
ro (:,2)= ro (:,1);
tau(:,2)= tau(:,1);
gdv(:,2)= gdv(:,1);
v  (:,2)= v  (:,1);
vrel(:,2)=vrel(:,1);
mur(:,2)= mur(:,1);
T  (:,2)= T  (:,1);



% 3.3 * Cálculo de las variables en el tiempo
%       -Se emplean derivadas temporales retrasadas y espaciales centradas.
for i=3:npt
%     ** Inicio del bucle temporal: Primera estimación de las variables
    [cp(:,i)]=calor_esp(cpini,roini,Tini,ro(:,i-1),T(:,i-1),p(i),k_cp1,...
        k_cp2,beta_cpini,b_cp1,b_cp2,npf);
    [k(:,i)]=conductividad(kini,T(:,i-1),p(i),C_k1,C_k2,npf);
    [mu(:,i)]=viscosidad_bv(muini,p(i),T(:,i-1),Tini,alfa_mu,beta_mu,npf);
    [eta(:,i)]=viscosidad_av(mu(:,i),gdv(:,i-1),G,npf,cCarr,expley);
    [ro(:,i)]=densidad(roini,p(i),T(:,i-1),Tini,epsilon(i),C_ro1,C_ro2,...
        npf);
    [tau(:,i)]=cortante(tau0,eta(:,i),gdv(:,i-1),npf);
    [gdv(:,i)]=gradvel(tau(:,i),tau0,eta(:,i),npf);
    [T(:,i)]=temperatura(T(:,i-1),Tinf,Tsup,ro(:,i),ro(:,i-1),cp(:,i),...
        ht,k(:,i),hz,tau(:,i),gdv(:,i),epsilon(i),dp(i),npf);
    
    
    % ** Cálculo iterativo de las variables para el instante i
    while max(difT)>error
        N(i)=N(i)+1;     % número de iteraciones en el instante i
        if N(i)==nit
            disp('superado nº max de iteraciones de gdv en paso:'),i
            pause    % Pulsar Ctrl+C en la command window para detenerlo
        else
            [cp(:,i)]=calor_esp(cpini,roini,Tini,ro(:,i),T(:,i),p(i),...
                k_cp1,k_cp2,beta_cpini,b_cp1,b_cp2,npf);
            [k(:,i)]=conductividad(kini,T(:,i),p(i),C_k1,C_k2,npf);
            [mu(:,i)]=viscosidad_bv(muini,p(i),T(:,i),Tini,alfa_mu,...
                beta_mu,npf);
            [eta(:,i)]=viscosidad_av(mu(:,i),gdv(:,i),G,npf,cCarr,expley);
            [ro(:,i)]=densidad(roini,p(i),T(:,i),Tini,epsilon(i),C_ro1,...
                C_ro2,npf);
            [tau(:,i)]=cortante(tau0,eta(:,i),gdv(:,i),npf);
            [gdv(:,i)]=gradvel(tau(:,i),tau0,eta(:,i),npf);
            [Tnew]=temperatura(T(:,i-1),Tinf,Tsup,ro(:,i),ro(:,i-1),...
                cp(:,i),ht,k(:,i),hz,tau(:,i),gdv(:,i),epsilon(i),dp(i),...
                npf);
            
            difT=abs(T(:,i)-Tnew)./T(:,i);   % diferencia entre iteraciones
            T(:,i)=Tnew;    % almacenamiento del valor calculado
        end
        [v(:,i)]=velocidad(gdv(:,i),hz,npf);
        vrel(:,i)=v(:,i)./v(:,1);
        [mur(:,i)]=tau(:,i)/p(i);
    end
    difT=(error+1)*ones(npf,1);   % puesta a cero del vector de diferencias
end
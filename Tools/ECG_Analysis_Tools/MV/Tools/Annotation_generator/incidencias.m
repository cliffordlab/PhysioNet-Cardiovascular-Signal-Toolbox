function [tm,ids,LOG,tn,cn]=incidencias(torg,Tol,limitacion)
% [tm,ids,LOG]=incidencias(torg,Tol)
% Detecta los latidos anormales (ectopic beat detection).
% tm  = Latidos normales
% tn  = Latidos con las correcciones EDU
% ids = Indices de los latidos posteriores a cada incidencia
% LOG = Informe en formato cadena de caracteres
% torg = instantes de ocurrencia de latidos (en seg)
% Tol = Tolerancia a ectopicos. Por defecto Tol=1. (Tolerancia a ectopicos. Por defecto Tol=1.)
% limitacion = indica si se hace analisis cuando hay variaciones bruscas HRV (0 hace analisis, 1 no hace analisis)
% (c) Javier Mateo, 7-Dic-2001.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

start=clock; RC=[char(13),char(10)];

if nargin<2, Tol=1;end
if nargin<3, limitacion=1;end
tn=torg(:);tn=[tn(diff(tn)>0);tn(end)];cn=repmat('N',length(tn),1); %En tn estan las detecciones QRS, elimina las que estan "desordenadas"
LOG = [  'Informe de anomalias en la serie de latidos',RC];
LOG=[LOG,'===========================================',RC,RC];

ddt=diff(1./diff(tn));
ddt2=ddt.^2; %Aplica un filtro de mediana ventana 180, sobre ddt2
wind=min(180,2*length(tn)-4);
mf=medfilt1([flipud(ddt2(1:wind/2));ddt2;flipud(ddt2(end-wind/2+1:end))],wind-1); %Añade un inicio y un fin " girados" para el calculo
%El valor del vector umbral es la raiz cuarta de ddt2
TOL=Tol*(mf(wind/2+1:end-wind/2)).^0.25;TOL(TOL<0.05)=0.05;TOL(TOL>0.5)=0.5;

if any(TOL>=0.5) & 0
    tm=[];ids=[];tn=[];cn=[];
    LOG=[LOG,'El fichero presenta zonas de extrema variabilidad.',RC];
    LOG=[LOG,'Se desaconseja el estudio de HRV en este fichero.',RC];
    return
end
%Mira que latido presentan incidencias
idx=find(abs(ddt)>TOL);
while ~isempty(idx) & idx(1)<3 %Si hay incidencias en los primeros latidos las elimina
    tn=tn(2:end);
    cn=cn(2:end);
    TOL=TOL(2:end);
    ddt=diff(1./diff(tn));
    idx=find(abs(ddt)>TOL);
end

while ~isempty(idx) %Analiza cada incidencia
    id=idx(1);
    N=1;
    T=TOL(id);
    C=0;
    while C==0 & length(tn)>id+N+2 %Mientras no se solvente y queden latidos va tomando el siguiente lat
        N=N+1;
        [C,tm]=check(tn(id-2:id+N+2),T);
    end
    if C
        tn=[tn(1:id-2+C);tm;tn(id+N-2+C:end)];
        cn=[cn(1:id-2+C);repmat('x',length(tm),1);cn(id+N-2+C:end)];
        TOL=[TOL(1:id-2+C);repmat(T,length(tm),1);TOL(id+N-2+C:end)];
    else %Cuando se termina el vector y no se ha resuelto => elimina latido
        tn=tn(1:end-1);
        cn=cn(1:end-1);
        TOL=TOL(1:length(tn)-2);
    end
    ddt=diff(1./diff(tn));
    idx=find(abs(ddt)>TOL);
end

%id=find(filtfilt([1 1],1,~ismember(cn,'Nb?'))<1.5);
id=find(cn=='N');
tm=tn(id); %pone en tm los latidos marcados como normales
ids=find(diff(id)>1)+1; %Marca los latidos posteriores a incidencias

procdur=etime(clock,start);

LOG=[LOG,'Posicion primer latido original: ',num2str(torg(1)),' s.',RC];
LOG=[LOG,'Posicion primer latido normal:   ',num2str(tm(1)),' s.',RC];
LOG=[LOG,'Posicion ultimo latido original: ',num2str(torg(end)),' s.',RC];
LOG=[LOG,'Posicion ultimo latido normal:   ',num2str(tm(end)),' s.',RC];
LOG=[LOG,'Numero de latidos originales:    ',num2str(length(torg)),RC];
LOG=[LOG,'Numero de latidos normales:      ',num2str(length(tm)),RC];
LOG=[LOG,'Numero de anomalias encontradas: ',num2str(length(ids)),RC];
LOG=[LOG,'Maximo intervalo entre normales: ',num2str(max(diff(tm))),' s.',RC];
LOG=[LOG,'Periodo cardiaco medio:          ',num2str(mean(diff(tn)),2),' s.',RC];
LOG=[LOG,'Duracion procesado de anomalias: ',num2str(procdur,2),' s.',RC];


function [C,tm]=check(tt,T)
% C indica el indice donde se debe insertar
% tm son los valores a insertar
%tt contiene los tiempos de ocurrencia de los latidos a analizar
for i=1:3
xi=tt(i+1); yi=xi-tt(i);
xf=tt(end-3+i); yf=xf-tt(end-4+i);
tm{i}=tgap(xi,yi,xf,yf);
tp=[tt(i:i+1);tm{i};tt(end-3+i)];
dd(i)=max(abs(diff(1./diff(tp))));
end
[dd,C]=min([T,dd]); %Compara si las modificaciones dd ahora cumplen el umbral
C=C-1;   %C=0 => T<dd NO RESUELTO   C>1 => T>dd Resuelto el problema, se cumple la condicion C tiene la informacion del bucle i(1,2o3) que la cumple
if C>0
    tm=tm{C}(1:end-1);
else
    tm=[];
end

function [xm,ym]=tgap(x0,y0,xf,yf)
m=(yf-y0)/(xf-x0); g=1/(1-m);
k=1;
b0=y0-m*x0;
b1=(xf-yf-x0*g^k)/sum(g.^(1:k));
while 1
   k=k+1;
   b2=(xf-yf-x0*g^k)/sum(g.^(1:k));
   db0=abs(b1-b0); db1=abs(b2-b0);
   if (db0-db1)<0, break, end
   b1=b2;
end
N=k-1;
xm(1)=(x0+b1)*g;
ym(1)=xm(1)-x0;
for i=2:N
   xm(i,1)=g*(xm(i-1)+b1);
   ym(i,1)=xm(i)-xm(i-1);
end

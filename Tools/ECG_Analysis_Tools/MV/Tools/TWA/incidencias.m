function [tm,ids,LOG,tn,cn]=incidencias(torg,Tol,limitacion)
% [tm,ids,LOG]=incidencias(torg,Tol)
% OVERVIEW, Detects abnormal beats (ectopic beat detection).
% tm  = Normal beats
% tn  = Beats with EDU corrections
% ids = Indices of the beats after each incident
% LOG = Report in character string format
% torg = moments of heartbeat occurrence (in sec)
% Tol = Ectopic tolerance. Default Tol=1.
% limitacion = indicates if analysis is done when there are abrupt variations in HRV (0 does analysis, 1 does not do analysis)
% (c) Javier Mateo, 7-Dic-2001.
%  
%   REPO:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:
%       Written by Javier Mateo
%	COPYRIGHT (C) Javier Mateo, 7-Dic-2001
%   LICENSE:
%       This software is offered freely and without warranty under
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.  


start=clock; RC=[char(13),char(10)];

if nargin<2, Tol=1;end
if nargin<3, limitacion=1;end
tn=torg(:);tn=[tn(diff(tn)>0);tn(end)];cn=repmat('N',length(tn),1); %In tn are QRS detections, remove those that are "messy"
LOG = [  'Informe de anomalias en la serie de latidos',RC];
LOG=[LOG,'===========================================',RC,RC];

ddt=diff(1./diff(tn));
ddt2=ddt.^2; %Apply a mediam window filter 180, on ddt2
wind=min(180,2*length(tn)-4);
mf=medfilt1([flipud(ddt2(1:wind/2));ddt2;flipud(ddt2(end-wind/2+1:end))],wind-1); %Add a "turned" start and end to the calculation
%The value of the threshold vector is the fourth root of ddt2
TOL=Tol*(mf(wind/2+1:end-wind/2)).^0.25;TOL(TOL<0.05)=0.05;TOL(TOL>0.5)=0.5;

if any(TOL>=0.5) & 0
    tm=[];ids=[];tn=[];cn=[];
    LOG=[LOG,'El fichero presenta zonas de extrema variabilidad.',RC];
    LOG=[LOG,'Se desaconseja el estudio de HRV en este fichero.',RC];
    return
end
%Look what beat have incidents
idx=find(abs(ddt)>TOL);
while ~isempty(idx) & idx(1)<3 %If there are incidents in the first beats it eliminates them
    tn=tn(2:end);
    cn=cn(2:end);
    TOL=TOL(2:end);
    ddt=diff(1./diff(tn));
    idx=find(abs(ddt)>TOL);
end

while ~isempty(idx) %Analyze each incident
    id=idx(1);
    N=1;
    T=TOL(id);
    C=0;
    while C==0 & length(tn)>id+N+2 %As long as it is not solvent and there are still beats, it takes the next lat
        N=N+1;
        [C,tm]=check(tn(id-2:id+N+2),T);
    end
    if C
        tn=[tn(1:id-2+C);tm;tn(id+N-2+C:end)];
        cn=[cn(1:id-2+C);repmat('x',length(tm),1);cn(id+N-2+C:end)];
        TOL=[TOL(1:id-2+C);repmat(T,length(tm),1);TOL(id+N-2+C:end)];
    else %When the vector is finished and has not been resolved => eliminates heartbeat
        tn=tn(1:end-1);
        cn=cn(1:end-1);
        TOL=TOL(1:length(tn)-2);
    end
    ddt=diff(1./diff(tn));
    idx=find(abs(ddt)>TOL);
end

%id=find(filtfilt([1 1],1,~ismember(cn,'Nb?'))<1.5);
id=find(cn=='N');
tm=tn(id); %tm beats marked as normal
ids=find(diff(id)>1)+1; %Mark the beats after incidents

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
% C indicates the index where it should be inserted
% tm are the values ​​to insert
% tt contains the times of occurrence of the beats to analyze
for i=1:3
xi=tt(i+1); yi=xi-tt(i);
xf=tt(end-3+i); yf=xf-tt(end-4+i);
tm{i}=tgap(xi,yi,xf,yf);
tp=[tt(i:i+1);tm{i};tt(end-3+i)];
dd(i)=max(abs(diff(1./diff(tp))));
end
[dd,C]=min([T,dd]); %Compare if dd modifications now meet the threshold
C=C-1;   %C = 0 => T <dd NOT SOLVED C> 1 => T> dd Once the problem is resolved, condition C has the information of loop i (1,2o3) that meets it
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

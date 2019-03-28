function heasig=readheader(name)

% function heasig=readheader(name);
%
% READHEADER function reads the header of DB signal files
%	Input parameter: character string with name of header file
%	Output parameter: struct heasig with header information
%
% Salvador Olmos
% e-mail: olmos@posta.unizar.es
% LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%

% Opening header file
fid=fopen(name,'rt');
if (fid<=0)
   disp(['error in opening file ' name]);
end

pp=' /+:()x';

% First line reading
s=fgetl(fid);
% Remove blank or commented lines
while s(1)=='#'
  s=fgetl(fid);
end

[heasig.recname,s]=strtok(s,pp);
[s1 s]=strtok(s,pp);
heasig.nsig=str2num(s1);
[s1 s]=strtok(s);
if isempty(findstr(s1,'/'))
   heasig.freq=str2num(s1);
else
   s1=strtok(s1,'/');
   heasig.freq=str2num(s1);
  
end

[s1 s]=strtok(s,pp);
heasig.nsamp=str2num(s1);

if ~isempty(deblank(s))
   [s1 s]=strtok(s,':');
   hour=str2num(s1);
   [s1 s]=strtok(s,':');
   min=str2num(s1);
   [s1 s]=strtok(s,pp);
   sec=str2num(s1);  
end

if ~isempty(deblank(s))
   [s1 s]=strtok(s,'/');
   day=str2num(s1);
   [s1 s]=strtok(s,'/');
   month=str2num(s1);
   [s1 s]=strtok(s,pp);
   year=str2num(s1);  
end
if exist('hour','var') heasig.date=datenum(year,month,day,hour,min,sec); end

% default values
for i=1:heasig.nsig
  heasig.units(i,:)='mV';
end
sig=1;

% Reading nsig lines, corresponding one for every lead
for i=1:heasig.nsig
  s=fgetl(fid);
  % Remove blank or commented lines
  while s(1)=='#'
    s=fgetl(fid);
  end

  [heasig.fname(i,:),s]=strtok(s,pp);
  [s1,s]=strtok(s,pp);
  if i==1  heasig.group(i)=0;
  else
     if strcmp(heasig.fname(i,:),heasig.fname(i-1,:)) 
        heasig.group(i)=0;
     else
        heasig.group(i)=heasig.group(i-1)+1;
     end
  end
  a=[findstr(s,'x') findstr(s,':') findstr(s,'+')];
  [s2,aux]=strtok(s); % rute prob com o bd do politecnico 13/2/03
  a(a>length(s2))=[]; % rute prob com o bd do politecnico 13/2/03
  if isempty(a)
     heasig.fmt(i)=str2num(s1);     
  else
    [s2,s]=strtok(s);
    a=[a length(s2)+1];     
    for k=1:length(a)-1
      switch (s2(a(k)))
       case '+',
     	heasig.fmt(i)=str2num(s1);
     	heasig.offset(i)=str2num(s2(a(k)+1:a(k+1)-1));  
       case ':',
     	heasig.fmt(i)=str2num(s1);
        heasig.skew(i)=str2num(s2(a(k)+1:a(k+1)-1));
       case 'x',
     	heasig.fmt(i)=str2num(s1);
     	heasig.spf(i)=str2num(s2(a(k)+1:a(k+1)-1));  
      end
    end
  end
  [s1,s]=strtok(s,pp);  
  a=[findstr(s,'(') findstr(s,'/')];
  if isempty(s1)
      heasig.gain(i)=0;
      heasig.baseline(i)=0;
  else
      if isempty(a)
        heasig.gain(i)=str2num(s1);
      else
       [s2,s]=strtok(s);
        %a=[a length(s2)+1];
        a=[a(1) length(s2)+1]; % rute prob com o bd do politecnico na unidade da respiraçao(dl/sec) 13/2/03 
        for k=1:length(a)-1
  	  switch (s2(a(k)))
           case '(', 
		heasig.gain(i)=str2num(s1);
		a2=findstr(s2,')');
		heasig.baseline(i)=str2num(s2(1+a(k):a2-1));
	   case '/',
		heasig.gain(i)=str2num(s1);
                f=s2(a(k)+1:end);
		heasig.units(i,1:length(f))=f;
          end
        end
      end
      [s1,s]=strtok(s,pp);
      heasig.adcres(i)=str2num(s1);
      [s1,s]=strtok(s,pp);
      heasig.adczero(i)=str2num(s1);
      [s1,s]=strtok(s,pp);
      heasig.initval(i)=str2num(s1);
      [s1,s]=strtok(s,pp);
      heasig.cksum(i)=str2num(s1);
      [s1,s]=strtok(s,pp);
      heasig.bsize(i)=str2num(s1);
      heasig.desc(i,1:length(s))=s;  
  end  
end

fclose(fid);

function  [anntime,type,subtype,chan,num,comments]=read_ann(recordName,annotator)
% Reads a WFDB annotation and returns:
%
% anntime
%       Nx1 vector of the ints. The time of the annotation in samples
%       with respect to the fist sample in the signals in WFDB dat file.
%       To convert this vector to a string of time stamps see WFDBTIME.
%
% type
%       Nx1 vector of the chars describing annotation type. 
%       For a list of standard type used by PhyioNet, please see:
%           https://www.physionet.org/physiotools/wfdb/lib/ecgcodes.h
%
% subtype
%       Nx1 vector of the ints describing annotation subtype, if not empty.
%       For a list of standard annotation codes used by PhyioNet, please see:
%             http://www.physionet.org/physiobank/annotations.shtml
%
% chan
%       Nx1 vector of the ints describing annotation subtype, if not empty.
%
% num
%       Nx1 vector of the ints describing annotation NUM, if not empty.
%
% comments
%       Nx1 vector of the cells describing annotation comments, if not empty.
%
%
% Required Parameters:
%
% recorName
%        String specifying the name of the record
%  
% annotator
%        String specifying the name of the annotation file 
%
%   ORIGINAL SOURCE AND AUTHORS:     
%       Qiao Li 
%       
%       Dependent scripts written by various authors 
%       (see functions for details)       
%	COPYRIGHT (C) 2016 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%
if isnumeric(recordName)
    recordName = num2str(recordName);
end

annfile=[recordName '.' annotator];
fid=fopen(annfile,'r');
A= fread(fid, [2, inf], 'uint8')';
fclose(fid);
anntime=[];
type=[];
chan=[];
num=[];
comments=[];
subtype_temp=[];
curnum=0;
curchn=0;
sa=size(A);
saa=sa(1);
i=1;

ann_n=0;
while i<=saa
    % Each annotation occupies an even number of bytes. The first byte in 
    % each pair is the least significant byte. The six most significant 
    % bits (A) of each byte pair are the annotation type code, and the ten 
    % remaining bits (I) specify the time of the annotation, measured in 
    % sample intervals from the previous annotation (or from the beginning 
    % of the record for the first annotation). If 0 < A <= ACMAX, then A 
    % is defined in <wfdb/ecgcodes.h>. (MIT format)
    typeh=bitshift(A(i,2),-2);
    if typeh==59 % SKIP, the next 4 bytes are the interval
        if bitshift(A(i+3,2),-2)~=0
            type=[type;int2ann(bitshift(A(i+3,2),-2))];
            anntime=[anntime;A(i+2,1)+bitshift(A(i+2,2),8)+...
                    bitshift(A(i+1,1),16)+bitshift(A(i+1,2),24)];
            ann_n=ann_n+1;
            num(ann_n)=curnum;
            chan(ann_n)=curchn;
        end
        i=i+3;
    elseif typeh==60 % NUM, I = annotation num field for current and subsequent annotations; otherwise, assume previous annotation num (initially 0). 
        % curnum=bitshift(bitand(A(i,2),3),8)+A(i,1);
        curnum=A(i,1);
        if curnum>127
            curnum=curnum-255-1;
        end
        num(ann_n)=curnum;
    elseif typeh==61 % SUB, I = annotation subtyp field for current annotation only; otherwise, assume subtyp = 0. 
        % cursub=bitshift(bitand(A(i,2),3),8)+A(i,1);
        cursub=A(i,1);
        if cursub>127
            cursub=cursub-255-1;
        end
        subtype_temp(ann_n)=cursub;
    elseif typeh==62 % CHN, I = annotation chan field for current and subsequent annotations; otherwise, assume previous chan (initially 0). 
        % curchn=bitshift(bitand(A(i,2),3),8)+A(i,1);
        curchn=A(i,1);
        if curchn>127
            curchn=curchn-255-1;
        end
        chan(ann_n)=curchn;
    elseif typeh==63 % AUX, I = number of bytes of auxiliary information (which is contained in the next I bytes); an extra null, not included in the byte count, is appended if I is odd. 
        hilfe=bitshift(bitand(A(i,2),3),8)+A(i,1);
        hilfe=hilfe+mod(hilfe,2);
        curaux=[];
        for j=1:hilfe/2
            curaux=[curaux char(A(i+j,1)) char(A(i+j,2))];
        end
        comments{ann_n}=curaux;
        i=i+hilfe/2;
    else
        anntime=[anntime;bitshift(bitand(A(i,2),3),8)+A(i,1)];
        type=[type;int2ann(bitshift(A(i,2),-2))];
        ann_n=ann_n+1;
        num(ann_n)=curnum;
        chan(ann_n)=curchn;
    end
   i=i+1;
end
type(length(type))=[];       % last line = EOF (=0)
anntime(length(anntime))=[];   % last line = EOF
num(length(num))=[];
chan(length(chan))=[];
num=num';
chan=chan';
if ~isempty(comments)
    comments=comments';
end
subtype=zeros(length(anntime),1);
if ~isempty(subtype_temp)
    subtype(find(subtype_temp~=0))=subtype_temp(find(subtype_temp~=0));
end
clear A;
anntime= (cumsum(anntime));

end

function ann_type=int2ann(ann_int)
% input: ann_int, annotation code, integer
% output:ann_type, annotation type, char

Typestr='NLRaVFJASEj/Q~|sT*D"=pB^t+u?![]en@xf(`)''r';
codeint=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,39,40,40,41];

findp=find(codeint==ann_int);
if ~isempty(findp)
    ann_type=Typestr(findp(1));
else
    ann_type=' ';
end
end


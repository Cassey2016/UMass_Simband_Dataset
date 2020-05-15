function [refECG_pkloc,wbwrefHR] = my_ECG_Peak_Detection_concise(output_data,fs,Ref_ECG_buffer)
% unload data
new_ECG_peak_loc = output_data.new_ECG_peak_loc;
my_Simband_datafilename = output_data.my_Simband_datafilename;
iiii_ECG_start = output_data.iiii_ECG_start;
iiii_ECG_end = output_data.iiii_ECG_end;
    
% subfunction 1
[refECG_pkloc,wbwrefHR] = my_ECG_ref_buffer(new_ECG_peak_loc,my_Simband_datafilename,Ref_ECG_buffer,fs,iiii_ECG_start,iiii_ECG_end);

end

% subfunction 1
function [refECG_pkloc,wbwrefHR] = my_ECG_ref_buffer(new_ECG_peak_loc,my_Simband_datafilename,Ref_ECG_buffer,fs,iiii_ECG_start,iiii_ECG_end)

    if ((strcmp(my_Simband_datafilename,'4005'))...
            || (strcmp(my_Simband_datafilename,'4002'))...
            || (strcmp(my_Simband_datafilename,'4001'))...
            || (strcmp(my_Simband_datafilename,'4006'))...
            || (strcmp(my_Simband_datafilename,'4012'))...
            || (strcmp(my_Simband_datafilename,'4016'))...
            || (strcmp(my_Simband_datafilename,'4013'))...
            || strcmp(my_Simband_datafilename,'4033')...
            || strcmp(my_Simband_datafilename,'4045')...
            || strcmp(my_Simband_datafilename,'4043')...
            || strcmp(my_Simband_datafilename,'4042')...
            || strcmp(my_Simband_datafilename,'4040')...
            || strcmp(my_Simband_datafilename,'4038')...
            || strcmp(my_Simband_datafilename,'4037')...
            || strcmp(my_Simband_datafilename,'4027')...
            || strcmp(my_Simband_datafilename,'4025')...
            || strcmp(my_Simband_datafilename,'4024')...
            || strcmp(my_Simband_datafilename,'4022')...
            || strcmp(my_Simband_datafilename,'4021')...
            || strcmp(my_Simband_datafilename,'4019')...
            || strcmp(my_Simband_datafilename,'4015'))

        temp_idx = find(new_ECG_peak_loc >= iiii_ECG_start & new_ECG_peak_loc <= iiii_ECG_end);
        refECG_pkloc = new_ECG_peak_loc(temp_idx) - iiii_ECG_start;
        wbwrefHR = 60./diff(refECG_pkloc)*fs;
    elseif (strcmp(my_Simband_datafilename,'Shiju_tightness02')...
        || strcmp(my_Simband_datafilename,'Shiju_tightness01'))
        temp_idx = find(new_ECG_peak_loc >= iiii_ECG_start & new_ECG_peak_loc <= iiii_ECG_end);
        refECG_pkloc = new_ECG_peak_loc(temp_idx) - iiii_ECG_start;
        wbwrefHR = 60./diff(refECG_pkloc)*fs;
    else
        [refECG_pkloc,wbwrefHR] = myECGpeakdetection(Ref_ECG_buffer,fs);
    end
    refECG_pkloc(refECG_pkloc < 1) = []; % change those possible zero index as one.
    refECG_pkloc(refECG_pkloc > 30*fs) = []; % remove those index greater than 30 second length.
end

function [refECG_pkloc,wbwrefHR] = myECGpeakdetection(ECG_input,fs)
    HDR = qrsdetect(ECG_input,fs,2);
    refECG_pkloc = HDR.EVENT.POS;
    wbwrefHR = 60*fs./(abs(diff(refECG_pkloc)));
end

function [H2,HDR,s] = qrsdetect(fn,arg2,arg3,varargin)
% QRSDETECT - detection of QRS-complexes
%
%   HDR = qrsdetect(fn,chan,Mode)
%   ... = qrsdetect(fn, 0, Mode, '-o',outputFilename)
%   ... = qrsdetect(... ,'-e',eventFilename)
%   HDR = qrsdetect(s,Fs,Mode) 
%
% INPUT
%   	fn	filename
%   	chan    channel number of ecg data
%		if chan is empty, all channels which contain ECG, ecg, or EKG in HDR.Label 
%		are used. 
%   	s       ecg signal data 
%   	Fs      sample rate 
%   	Mode    optional - default is 2
%               1: method [1] is used
%		2: method [2] is used
%	outputFilename
%		name of file for storing the resulting data with the
%		detected spikes and bursts in GDF format.
%		
%	eventFilename
%		filename to store event inforamation in GDF format. this is similar to 
%		the outputFile, except that the signal data is not included and is, therefore,
%		much smaller than the outputFile
%
% OUTPUT
%   HDR.EVENT  fiducial points of qrs complexes	
%
%
% see also: PROCESSING, EVENTCODES.TXT, SLOAD 
%
% Reference(s):
% [1] M.-E. Nygards, L. Sörnmo, Delineation of the QRS complex using the envelope of the e.c.g
%       Med. & Biol. Eng. & Comput., 1983, 21, 538-547.
% [2] V. Afonso, W. Tompkins, T. Nguyen, and S. Luo, "ECG beat detection using filter banks."
% 	IEEE Trans. Biomed. Eng. 46(2):192-202, Feb. 1999.

% $Id$
% Copyright (C) 2000-2003,2006,2009,2011 by Alois Schloegl <alois.schloegl@gmail.com>
% This is part of the BIOSIG-toolbox http://biosig.sf.net/

% BIOSIG is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 3 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


chan = 0; 
MODE = 2;       % default
if (nargin==2) 
	if isnumeric(arg2)
		chan = arg2;
	else
		MODE = arg3;
	end;
elseif (nargin>=3) 
	chan = arg2;
	MODE = arg3;
end;

%%%%% default settings %%%%%
outFile = [];
evtFile = [];

%%%%% analyze input arguments %%%%%
k = 1;
while k <= length(varargin)
	if ischar(varargin{k})
		if (strcmp(varargin{k},'-o'))
			k = k + 1;
			outFile = varargin{k};
		elseif (strcmp(varargin{k},'-e'))
			k = k + 1;
			evtFile = varargin{k};
		else
			warning(sprintf('unknown input argument <%s>- ignored',varargin{k}));
		end;
	end;
	k = k+1;
end;


if isnumeric(fn),
	s = fn;
	HDR.SampleRate = arg2;
	HDR.NS = size(s,2);
	chan = 1:size(s,2);
	CHAN = 1:size(s,2);
else
	if ~isempty(outFile) 
		[s,HDR] = sload(fn);
		if chan==0 || isempty(chan);
			chan = unique([strmatch('ecg',HDR.Label),strmatch('ECG',HDR.Label),strmatch('EKG',HDR.Label)]);
		end;	
		CHAN = chan;
	else 
		HDR = sopen(fn,'r',chan);HDR = sclose(HDR);
		if chan==0 || isempty(chan);
			chan = unique([strmatch('ecg',HDR.Label),strmatch('ECG',HDR.Label),strmatch('EKG',HDR.Label)]);
		end;	
		[s,HDR] = sload(fn,chan);
		CHAN = chan;
		chan = 1:length(CHAN);
	end; 
end;

ET = [];
for ch = 1:length(chan)',
	k = chan(ch);

	if 0,

        elseif MODE==2,     % QRS detection based on Afonso et al. 1999
                POS = nqrsdetect(s(:,k),HDR.SampleRate);

        elseif MODE==1,     % QRS detection based on Hilbert transformation. For details see [1]
                y   = processing({'ECG_envelope',HDR.SampleRate},s(:,k));
                TH  = quantile(y,.90);

                POS = gettrigger(y,TH);	% fiducial point

        else %if Mode== ???     % other QRS detection algorithms
		fprintf(2,'Error QRSDETECT: Mode=%i not supported',Mode);
		return;
        end;

        % difference between R-peak and fiducial point
        [t,sz] = trigg(s(:,k),POS,floor(-HDR.SampleRate),ceil(HDR.SampleRate));
        [tmp,ix] = max(abs(mean(reshape(t,sz(2:3)),2)));
        delay = HDR.SampleRate + 1 - ix;

        ET = [ET; [POS-delay, repmat([hex2dec('0501'),CHAN(ch),0], size(POS,1),1)]];
end;


[tmp,ix] = sort(ET(:,1));
if isfield(HDR,'T0')
	H2.T0 = HDR.T0;
end;
if isfield(HDR,'Patient')
	H2.Patient = HDR.Patient;
end;
H2.EVENT.POS = ET(ix,1);
H2.EVENT.TYP = ET(ix,2);
H2.EVENT.CHN = ET(ix,3);
H2.EVENT.DUR = ET(ix,4);
H2.EVENT.SampleRate = HDR.SampleRate;
H2.TYPE = 'EVENT';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(outFile)
		%%% write data to output
		HDR.TYPE  = 'GDF';
		HDR.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		HDR.FILE = [];
		HDR.FileName  = outFile;
		HDR.FILE.Path = '';
		HDR.PhysMax   = max(s);
		HDR.PhysMin   = min(s);
		HDR.DigMax(:) =  2^15-1;
		HDR.DigMin(:) = -2^15;
		HDR.GDFTYP(:) = 3;
		HDR.FLAG.UCAL = 0;
		HDR.NRec = size(s,1);
		HDR.SPR = 1;
		HDR.NS  = size(s,2);
		HDR.Dur = 1/HDR.SampleRate;
		HDR = rmfield(HDR,'AS');
		HDR.EVENT = H2.EVENT; 
		HDR = sopen(HDR,'w');
		if (HDR.FILE.FID < 0) 
			fprintf(2,'Warning can not open file <%s> - GDF file can not be written\n',HDR.FileName);
		else
			HDR = swrite(HDR,s);
			HDR = sclose(HDR);
		end; 
end;

if ~isempty(evtFile)
		%%% write data to output
		H2.TYPE  = 'EVENT';
		H2.VERSION = 3;
		%[p,f,e]=fileparts(fn);
		H2.FILE = [];
		H2.FileName  = evtFile;
		H2.FILE.Path = '';
		H2.NRec = 0;
		H2.SPR = 0;
		H2.Dur = 1/HDR.SampleRate;
		%H2 = rmfield(HDR,'AS');
		H2 = sopen(H2, 'w');
		if (H2.FILE.FID<0) 
			fprintf(2,'Warning can not open file <%s> - EVT file can not be written\n',H2.FileName); 
		else
			H2 = sclose(H2);
 		end;
end;


end

function QRS = nqrsdetect(S,fs);

% nqrsdetect - detection of QRS-complexes
%
%   QRS=nqrsdetect(S,fs);
%
% INPUT
%   S       ecg signal data
%   fs      sample rate
%
% OUTPUT
%   QRS     fiducial points of qrs complexes
%
%
% see also: QRSDETECT
%
% Reference(s):
% [1]: V. Afonso, W. Tompkins, T. Nguyen, and S. Luo, "ECG beat detection using filter banks,"
% 	IEEE Trans. Biomed. Eng., vol. 46, no. 2, pp. 192--202, Feb. 1999
% [2]: A.V. Oppenheim, R.W. Schafer, and J.R. Buck,  Discrete-Time Signal
% 	Processing, second edition, Prentice Hall, 1999, chapter 4.7.3

% Copyright (C) 2006 by Rupert Ortner, Graz University of Technology, Austria.
% Copyright (C) 2009 by Jens Stampe Soerensen, Aalborg University, Denmark.
%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, ...
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
%% USA

S=S(:);
S=full(S);
N=round(fs);   %Filter orde
%---------------------------------------
%Replaces filter bank in [1]
Bw=5.6;     %filter bandwidth
Bwn=1/(fs/2)*Bw;    
M=round((fs/2)/Bw); %downsampling rate

Wn1=[Bwn 2*Bwn];    %bandwidth of the second filter
Wn2=[2*Bwn 3*Bwn];
Wn3=[3*Bwn 4*Bwn];
Wn4=[4*Bwn 5*Bwn];

h1=fir1(N,Wn1,'DC-0');
h2=fir1(N,Wn2,'DC-0');
h3=fir1(N,Wn3,'DC-0');
h4=fir1(N,Wn4,'DC-0');

%Polyphase implementation of the filters
y=cell(1,5);    
y{2}=polyphase_imp(S,h1,M); %W1 
y{3}=polyphase_imp(S,h2,M); %W2
y{4}=polyphase_imp(S,h3,M); %W3
y{5}=polyphase_imp(S,h4,M); %W4
%----------------------------------------------


cut=ceil(N/M);  %Cutting off of initial transient because of the filtering
y2=[zeros(cut,1);y{2}(cut:length(y{2}))];
y3=[zeros(cut,1);y{3}(cut:length(y{3}))];
y4=[zeros(cut,1);y{4}(cut:length(y{4}))];
y5=[zeros(cut,1);y{5}(cut:length(y{5}))];
%----------------------------------------


P1=sum([abs(y2) abs(y3) abs(y4)],2); %see [1] equation (13)
P2=sum([abs(y2) abs(y3) abs(y4) abs(y5)],2);
P4=sum([abs(y3) abs(y4) abs(y5)],2);

FL1=MWI(P1); %Feature 1 according to Level 1 in [1]
FL2=MWI(P2); %Feature 2 according to Level 2
FL4=MWI(P4); %Feature 4 according to Level 4


%--------------------------------------
%Level 1 [1]
d=sign(diff(FL1));
d1=[0;d];
d2=[d;0];
f1=find(d1==1);
f2=find(d2==-1);
EventsL1=intersect(f1,f2); %Detected events

%-------------------------------------------------------
%Level 2 [1]
meanL1=sum(FL2(EventsL1),1)/length(EventsL1);
NL=meanL1-meanL1*0.1;   %Start Noise Level
SL=meanL1+meanL1*0.1;   %Start Signal Level
threshold1=0.08;    %Threshold detection block 1
threshold2=0.7;     %Threshold detection block 2
[DS1,Class1]=detectionblock(FL2,EventsL1,NL,SL,threshold1);
[DS2,Class2]=detectionblock(FL2,EventsL1,NL,SL,threshold2);

%---------------------------------------------------
%Level 3 [1]
ClassL3 = zeros(length(EventsL1),1);
for i=1:length(EventsL1)
    C1=Class1(i);
    C2=Class2(i);
    if C1==1
        if C2==1
            ClassL3(i) = 1;
        else
            delta1=(DS1(i)-threshold1)/(1-threshold1);
            delta2=(threshold2-DS2(i))/threshold2;
            if delta1>delta2
                ClassL3(i) = 1;
			else
				ClassL3(i) = 0;
            end
        end
    else
        if C2==1;
            ClassL3(i) = 1;
		else
			ClassL3(i) = 0;
        end
    end
end
SignalL3=EventsL1(ClassL3==1);   %Signal Level 3
NoiseL3=EventsL1(ClassL3==0); %Noise Level 3

%--------------------------------------------
%Level 4 [1]
threshold=0.3;
SL=(sum(FL4(SignalL3),1))/length(SignalL3);   %Initial Signal Level
NL=(sum(FL4(NoiseL3),1))/length(NoiseL3);      %Initial Noise Level

SignalL4 = zeros(1,length(SignalL3));
NoiseL4 = zeros(1,length(NoiseL3));
DsL4 = zeros(length(EventsL1),1);
for i=1:length(EventsL1)
    Pkt=EventsL1(i);    
    if ClassL3(i)==1;   %Classification after Level 3 as Signal
	   SignalL4(i) = Pkt;
	   SL=history(SL,FL4(Pkt));       
       Ds=(FL4(Pkt)-NL)/(SL-NL);       %Detection strength
       if Ds<0
           Ds=0;
       elseif Ds>1
           Ds=1;
	   end
	   DsL4(i) = Ds;
    else        %Classification after Level 3 as Noise
       Ds=(FL4(Pkt)-NL)/(SL-NL);   
       if Ds<0
           Ds=0;
       elseif Ds>1
           Ds=1;
	   end
	   DsL4(i) = Ds;       
       if Ds>threshold          %new classification as Signal  
           SignalL4(i) = Pkt;
		   SL=history(SL,FL4(Pkt));     
       else                      %new classification as Noise
           NoiseL4(i) = Pkt;
		   NL=history(NL,FL4(Pkt));      
       end
   end
end

SignalL4=SignalL4(SignalL4~=0);
NoiseL4=NoiseL4(NoiseL4~=0);
%Clean up;

%------------------------------------------------
%Level 5  
%if the time between two RR complexes is too long => go back and check the
%events again with lower threshold
SignalL5=SignalL4;
periods=diff(SignalL4);
M1=100;
a=1;
b=1/(M1)*ones(M1,1);
meanperiod=1.5*filter(b,a,periods); %mean of the RR intervals
SL=sum(FL4(SignalL4))/length(SignalL4);
NL=sum(FL4(NoiseL4))/length(NoiseL4);
threshold=0.2;
for i=1:length(periods)
    if periods(i)>meanperiod   %if RR-interval is to long
        intervall=SignalL4(i):SignalL4(i+1);
        critical=intersect(intervall,NoiseL4);   
        for j=1:length(critical)
            Ds=(FL4(critical(j))-NL)/(SL-NL); 
            if Ds>threshold         %Classification as Signal
                SignalL5=union(SignalL5,critical(j));				 
            end
        end
    end
end

clear meanperiod;
%---------------------------------------------------
%Umrechnung auf Originalsignal (nicht downgesamplet)
Signaln=conversion(S,FL2,SignalL5,M,N,fs);
%----------------------------------------------------
%Level 6 If interval of two RR-complexes <0.24 => go back and delete one of them
height=FL2(SignalL5);   
Signal=Signaln;
temp=round(0.1*fs);
difference=diff(Signaln);  %Difference between two signal points
k=find(difference<temp);
for i=1:length(k)
    verg=[height(k(i)),height(k(i)+1)]; 
    [x,j]=max(verg);    
    if j==1
        Signal=setxor(Signal,Signaln(k(i)+1)); %Deleting first Event
    else
        Signal=setxor(Signal,Signaln(k(i))); %Deleting second Event
    end
end
QRS=Signal(:);
end

%-------------------------------------------------------------------
%-------------------------------------------------------------------
%-------------------------------------------------------------------
%subfunctions

function y=MWI(S)

% MWI - Moving window integrator, computes the mean of two samples
%   y=MWI(S)
%
% INPUT
%   S       Signal
%
% OUTPUT
%   y       output signal
a=[0;S];
b=[S;0];
c=[a,b];
y=sum(c,2)/2;
y=y(1:length(y)-1);

end
%------------------------------------------------
function y=polyphase_imp(S,h,M)

% polyphase_imp - polyphase implementation of decimation filters [2]
%   y=polyphase_imp(S,h,M)
%
% INPUT
%   S       ecg signal data
%   h       filter coefficients
%   M       downsampling rate
%
% OUTPUT
%   y       filtered signal
%

%Determining polyphase components ek
e=cell(M,1);
l=1;
m=mod(length(h),M);
while m>0
    for n=1:ceil(length(h)/M)
        el(n)=h(M*(n-1)+l);
    end
    e{l}=el;
    l=l+1;
    m=m-1;
end
clear el;
for i=l:M
    for n=1:floor(length(h)/M)
        el(n)=h(M*(n-1)+i);      
    end
    e{i}=el;
end
%Filtering
max=ceil((length(S)+M)/M); 
w = zeros(max,M);
Sdelay = zeros(M-1+length(S),1);
Sdelay(M:end) = S; 
for i=1:M
	Sd= Sdelay(M-i+1:M:end);
	a=filter(e{i},1,Sd);
    if length(a)<max
        a=[a;zeros(max-length(a),1)]; 
    end
    w(:,i)=a;	
end
y=sum(w,2);
end
%----------------------------------------------------------
function [VDs,Class]=detectionblock(mwi,Events,NL,SL,threshold)


% detectionblock - computation of one detection block 
%
%   [Signal,Noise,VDs,Class]=detectionblock(mwi,Events,NL,SL,threshold)
%
% INPUT
%   mwi         Output of the MWI
%   Events      Events of Level 1 (see [1])
%   NL          Initial Noise Level
%   SL          Initial Signal Level
%   threshold   Detection threshold (between [0,1])
%
% OUTPUT
%   Signal      Events which are computed as Signal
%   Noise       Events which are computed as Noise
%   VDs         Detection strength of the Events
%   Class       Classification: 0=noise, 1=signal

signalCounter = 0;
noiseCounter = 0;
VDs=zeros(length(Events),1);
Class=zeros(length(Events),1);
sumsignal=SL;
sumnoise=NL;
for i=1:length(Events)
    P=Events(i);
    Ds=(mwi(P)-NL)/(SL-NL); %Detection strength	    
	if(Ds<0)
        Ds=0;    	
	elseif Ds>1
        Ds=1;
	end
    VDs(i) = Ds;
	
	if(Ds>threshold)     %Classification as Signal
		signalCounter = signalCounter + 1;
        Class(i) = 1;
		sumsignal=sumsignal+mwi(P);
		SL=sumsignal/(signalCounter+1);    %Updating the Signal Level
    else        %Classification as Noise
		noiseCounter = noiseCounter + 1;
        Class(i) = 0;
		sumnoise=sumnoise+mwi(P);
        NL = sumnoise/(noiseCounter + 1);		
	end	
end

end
%------------------------------------------------------------
function [pnew]=conversion(S,FL2,Signaln,M,N,fs)

% conversion - sets the fiducial points of the downsampled Signal on the
% samplepoints of the original Signal
% 
%   [pnew]=conversion(S,FL2,pold,M,N,fs)
%
% INPUT
%   S           Original ECG Signal
%   FL2         Feature of Level 2 [1]
%   pold        old fiducial points
%   M           M downsampling rate
%   N           filter order
%   fs          sample rate
%
% OUTPUT
%   pnew        new fiducial points
%

P=M;
Q=1;
FL2res=resample(FL2,P,Q);       %Resampling
%nans1=isnan(S);
%nans=find(nans1==1);
S(isnan(S)==1)=mean(S);    %Replaces NaNs in Signal
Signaln1 = Signaln + (M-1).*(Signaln-1);
%------------------- Sets the fiducial points on the maximum of FL2
Signaln2=Signaln1;  
Signaln2=Signaln2';     
int=2*M;    %Window length for the new fiducial point
Signaln3 = zeros(length(Signaln2),1);
for i=1:length(Signaln2)
    start=Signaln2(i)-int/2;
    if start<1
        start=1;
    end
    stop=Signaln2(i)+int/2;
    if stop>length(FL2res)
        stop=length(FL2res);
    end
    intervall=start:stop;       %interval
    FL2int=FL2res(intervall);
    pkt=find(FL2int==max(FL2int));  %Setting point on maximum of FL2
    if isempty(pkt)
		pkt=Signaln2(i)-start;
    else
        pkt=pkt(1); 
    end
    delay=N/2+M;
    Signaln3(i)=pkt+Signaln2(i)-int/2-delay;    %fiducial points according to FL2
end
%Sets the fiducial points on the maximum or minimum
%of the signal
Bw=5.6;   
Bwn=1/(fs/2)*Bw;
Wn=[Bwn 5*Bwn];
N1=32;
b=fir1(N1,Wn,'DC-0');
Sf=filtfilt(b,1,S);     %Filtered Signal with bandwidth 5.6-28 Hz
beg=round(1.5*M);
fin=1*M;

Signaln4 = zeros(length(Signaln3),1);
for i=1:length(Signaln3)
    start=Signaln3(i)-beg; % HD comment this. 07/20/2017
%     % HD add this part 07/20/2017
%     start = floor(Signaln3(i))-beg;
%     %
    if start<1
        start=1;
    end
    stop=Signaln3(i)+fin;
    if stop>length(Sf)
        stop=length(Sf);
    end
    intervall=start:stop;   %Window for the new fiducial point
	Sfint = abs(Sf(intervall)-sum(Sf(intervall))/length(Sf(intervall)));
    pkt=find(Sfint==max(Sfint));    %Setting point on maximum of Sfint
    if isempty(pkt)
        pkt=Signaln3(i)-start;
    else
        pkt=pkt(1); 
    end
    pkt=pkt(1);
    Signaln4(i)=pkt+Signaln3(i)-beg-1;
end
Signal=Signaln4';   %New fiducial points according to the original signal

fpointsb=Signal(Signal<N);
fpointse=Signal(Signal>length(S)-N);
pnew=setxor(Signal,[fpointsb fpointse]);

end
%-------------------------------------------
function yn=history(ynm1,xn)

% history - computes y[n]=(1-lambda)*x[n]+lambda*y[n-1]
%
%   yn=history(ynm1,xn)

lambda=0.8; %forgetting factor
yn=(1-lambda)*xn+lambda*ynm1;
end
function [x,sz] = trigg(s,TRIG,pre,post,gap)
% TRIGG cuts continous sequence into segments.
% Missing values (in case s is to short) are substituted by NaN's.  
%
% [X,sz] = trigg(s, TRIG, PRE, PST [, GAP])
%
% INPUT:
%  S	is the continous data sequence (1 channel per column)
%  TRIG	defines the trigger points
%  PRE 	offset of the start of each segment (relative to trigger) 
%  PST 	offset of the end of each segment (relative to trigger) 
%  GAP	number of NaN's to separate trials (default=0)
%  	TRIG, pre, post and gap are counted in samples
%
% OUTPUT:
%  X	is a matrix of size [sz(1), sz(2)*sz(3)]
%  sz	[size(s,2), post-pre+1+gap, length(TRIG)]   
%	sz(1) is the number of channels NS 
%	sz(2) is the number of samples per trial 
%	sz(3) is the number of trials i.e. length(TRIG)
%
% X3D = reshape(X,sz) returns a 3-dimensional matrix 
% X2D = reshape(X(K,:),sz(2:3)) returns channel K in a 2-dimensional matrix 
%
% see also: GETTRIGGER

% 	$Id$
%	Copyright (c) 1999-2005 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if nargin<5,
	gap = 0;
end;

post=round(post);

[nr,nc] = size(s);

% include leading nan's
off  = min(min([TRIG(:);+Inf])+pre-1,0);
% include following nan's
off2 = max(max([TRIG(:);-Inf])+post-length(s),0);
if ((off~=0) || (off2~=0))
	s    = [repmat(nan,-off,nc);s;repmat(nan,off2,nc)];
	TRIG = TRIG-off;
end; 

% devide into segments
N   = post-pre+1+gap;
sz  = [nc, post-pre+1+gap, length(TRIG)];
x   = repmat(NaN, [sz(1), sz(2)*sz(3)]);   
for m = 1:length(TRIG),
	x(:,m*N + (1-N:-gap)) = s(TRIG(m)+(pre:post)',:).';
end;

end
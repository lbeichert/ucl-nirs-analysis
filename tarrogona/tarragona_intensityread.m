function [result,stindex,endindex]=tarragona_intensityread(filename,average,stindex,endindex)
%Function to read intensity data from a Tarragona data file given by
% argument 'filename'. 
%
%   result=tarragona_intensityread(filename,average) 
% 
% The intensity data is averaged (as log intensity) using the argument
% 'average' as the number over which to average, (enter 0 for no averaging)
% and the data is returned as intensity and attenuation in OD.
%
% Created: 28/04/02
% Modified: 17/08/02 (comments only)
% (Veronica Hollis) 

% modified by Terence Leung last updated 26/5/2004
% stindex: starting index (sample) e.g. stindex=7 --> read from the 7
% spectrum
% endindex: ending index (samples)




if average<0
    error('Error: cannot have a negative value for "average" input argument');
end


% added Terence Leung
%------------------------------------
if nargin<4,
%     fprintf('\nTotal no. of samples : %i',spectra);
     stindex= input('\nstart sample :');
     endindex= input('\nend sample :');
end
% strpt= stindex;
% endpt= endindex;
Nspectra=endindex-stindex+1;
%---------------------------------



fid=fopen(filename,'r','ieee-le');      % open file for reading
fseek(fid,32,0);        % skip first record, no relevant information

%result.sd1version=fread(fid,1,'int32');
%result.sd1numfit=fread(fid,1,'int32');
%result.sd1extra=fread(fid,24,'uint8');

fseek(fid,10,0);        % skip first object (name) in next record

%result.sd2Nameindx=fread(fid,1);
%sd2Name=fread(fid,result.sd2Nameindx,'schar');
%result.sd2Name=char(sd2Name');
%fseek(fid,9-result.sd2Nameindx,0);

version=fread(fid,1,'uint16');      % program version, for future reference
datapoints=fread(fid,1,'uint16');        % read in number of data points
if average>datapoints
    error('"Average" argument exceeds number of spectra in file');
end
fseek(fid,4,0);     % skip next 4 bytes

%result.sd2ne=fread(fid,1,'int16');
%result.sd2np=fread(fid,1,'int16');

paths=fread(fid,1,'int16');     % number of pathlengths used
maths=fread(fid,1,'int16');     % number of mathematical operations
attn=fread(fid,1,'int16');      % number of attenuation
channels=fread(fid,1,'int16');      % number of channels
fseek(fid,6,0);     % skip info in next 6 bytes

%result.sd2hassd=fread(fid,1);
%result.sd2ccdtype=fread(fid,1,'uint8');
%result.sd2black=fread(fid,1,'int32');

extn=fread(fid,1,'int16');      % number of externals
fseek(fid,22,0);        % skip

%result.sd2sampleindx=fread(fid,1);
%sd2sample=fread(fid,result.sd2sampleindx,'schar');
%result.sd2sample=char(sd2sample');
%fseek(fid,21-result.sd2sampleindx,0);

datablock=fread(fid,8,'uint8')';        % array showing which CCD data blocks are used
fseek(fid,4,0);     % skip

%result.sd2xchanneldata=fread(fid,1,'uint8');
%result.sd2xchannelmode=fread(fid,1,'uint8');
%result.sd2pad1=fread(fid,2,'uint8');

events=fread(fid,1,'int32');        % number of events
fseek(fid,184,0);       % skip

%result.sd2Extra=fread(fid,184,'uint8')';

%--------------------------------------------------------------
fseek(fid,(channels+1)*345,0);      % skip channel data
    


% result.cd1Nameindx=fread(fid,1);
% cd1Name=fread(fid,result.cd1Nameindx,'schar');
% result.cd1Name=char(cd1Name')
% fseek(fid,255-result.cd1Nameindx,0);
% result.cd1colour=fread(fid,4);
% result.cd1used=fread(fid,1);
% result.cd1hassd=fread(fid,1);
% result.cd1linestyle=fread(fid,1);
% result.cd1linewidth=fread(fid,1);
% result.cd1pointstyle=fread(fid,1);
% result.cd1extra=fread(fid,15);
% result.cd1notview=fread(fid,1);
% result.cd1measurementindx=fread(fid,1);
% cd1measurement=fread(fid,result.cd1measurementindx,'schar');
% result.cd1measurement=char(cd1measurement');
% fseek(fid,31-result.cd1measurementindx,0);
% result.cd1unitindx=fread(fid,1);
% cd1unit=fread(fid,result.cd1unitindx,'schar');
% result.cd1unit=char(cd1unit');
% fseek(fid,31-result.cd1unitindx,0);

%------------------------------------------------------
if paths>0
    fseek(fid,paths*288,0);     % skip pathlength data if it exists
end    
    
%result.pdNameindx=fread(fid,1);
%pdName=fread(fid,result.pdNameindx,'schar');
%result.pdName=char(pdName');
%fseek(fid,255-result.pdNameindx,0);
%result.pdChannel=fread(fid,1,'int16');
%result.pdType=fread(fid,1,'int16');
%result.pdFP=fread(fid,1,'int16');
%result.pdLP=fread(fid,1,'int16');
%result.pdValue=fread(fid,1,'float');
%result.pdExtra=fread(fid,20);

if maths>0
    fseek(fid,maths*30,0);      % skip maths data if it exists
end

%result.md1channel=fread(fid,1,'int16');
%result.md1srceA=fread(fid,1,'int16');
%result.md1srceB=fread(fid,1,'int16');
%result.md1math=fread(fid,1,'int16');
%result.md1norm=fread(fid,1,'int16');
%result.md1pad=fread(fid,1,'int16');
%result.md1A0=fread(fid,1,'float');
%result.md1A1=fread(fid,1,'float');
%result.md1Extra=fread(fid,10,'uint8');

if attn>0
    fseek(fid,attn*28,0);     % skip attn data if it exists
end

if extn>0
    fseek(fid,extn*32,0);       % skip extn data if it exists
end

if events>0
%    fseek(fid,events*512); % skip events data if it exists(original)
    fseek(fid,events*512,0);  % skip events data if it exists (corrected 8/11/02)

end


% modified by Terence Leung (6/9/2004)
%------------------------------------------------------------------

% result.chan(:,8)  --> synch
%       .chan(:,9)  --> SaO2
%       .chan(:,10) --> POwave
%       .chan(:,11) --> portapres
%       .chan(:,14) --> syn

result.chan=zeros(Nspectra,64);
fseek(fid,512*(stindex+1),0);                          % move to the starting index
for i=1:Nspectra,
    result.chan(i,:)=fread(fid,64,'float')';           % contains channel data e.g. SaO2, POwave, synch                                                          
    fseek(fid,512-64*4,0);
end

fseek(fid,(datapoints-endindex-1)*512,0);
ftell(fid)
% fseek(fid,(datapoints+1)*512,0);        % skip info on each data point
% ftell(fid)
%result.dp1Channel=fread(fid,64,'float')';
%result.dp1Present=fread(fid,64)';
%result.dp1Datablock=fread(fid,8,'int32')';
%result.dp1Spare=fread(fid,5,'int32')';
%result.dp1Event=fread(fid,1,'int32')';
%result.dp1Triggertime=fread(fid,8,'int32')';
%result.dp1Trigger=fread(fid,8)';
%result.dp1extra=fread(fid,96)';

%-----------------------------------------------------

tic,
for n=1:8
    if datablock(n)==1
        Block(n).number=n;
        fseek(fid,4,0);     % skip
        %result.sd3Blacklevel=fread(fid,1,'int32');
        Block(n).Region=fread(fid,1,'int32');
        fseek(fid,23,0);
        %result.sd3Spare=fread(fid,23,'uint8');
        dataformat=fread(fid,1);
        if dataformat==1
            format='int32';
        end
        if dataformat==2
            format='float';
        end
        
        if dataformat==0        % added (28/1/2003 Terence)
            format='int32';
        end
        
        fseek(fid,128,0);
        %result.sh1Nameindx=fread(fid,1);
        %sh1Name=fread(fid,result.sh1Nameindx,'schar');
        %result.sh1Name=char(sh1Name');
        %fseek(fid,127-result.sh1Nameindx,0);
        suffixindx=fread(fid,1);
        suffix=fread(fid,suffixindx,'schar');
        Block(n).suffix=char(suffix');
        fseek(fid,127-suffixindx,0);
        pixels=fread(fid,1,'int32');
        spectra=fread(fid,1,'int32')
        fseek(fid,20,0);
        %result.sh1DataSze=fread(fid,1,'int32');
        %result.sh1Version=fread(fid,1,'int32');
        %result.sh1Extra=fread(fid,12);
        Block(n).wavelength=fread(fid,pixels,'float');
       Block(n).reference=fread(fid,pixels,format);
        
        fseek(fid,4*pixels,0);
        
% Re-marked by Terence 29/1/2003        
%         for i=1:spectra
%             Block(n).spectra(:,i)=fread(fid,pixels,format);
%         end
        
        % added Terence 26/5/04
%-------------------------------------------     

        fseek(fid,(stindex-1)*pixels*4,0);
        Block(n).spectra=fread(fid,[pixels,Nspectra],format);     % Terence 29/1/2003
%---------------------------------------------

        fseek(fid,spectra*2,0);
    end
end
% , toc

fclose(fid);            % close file

result.filename=filename;
result.num_spectra=spectra;
result.average=average;
result.wavelength=Block(1).wavelength;

% tic,
a=size(Block,2);
if average==0 || average==1
   for n=1:a
      %result.Block(n).reference=10.^Block(n).reference;
      %result.Block(n).intensity=10.^Block(n).spectra;

      % added 25/5/04 TSL
      result.Block(n).intensity=Block(n).spectra;
   end
      % commented 25/5/04 TSL
      % result.Block(n).attenuation=repmat(Block(n).reference(:),1,spectra) - Block(n).spectra;
else
   b=floor(spectra./average); 

   for n=1:a
     
      result.Block(n).reference=10.^Block(n).reference;
      result.Block(n).intensity=zeros(pixels,b);   % added by Terence 14/2/2003
      result.Block(n).attenuation=zeros(pixels,b);   % added by Terence 14/2/2003

      for j=1:b
         % changed by Terence Leung 13/10/2004
         result.Block(n).intensity(:,j)=((1./average).*sum(Block(n).spectra(:,j*average -average +1:j*average),2));

%          result.Block(n).intensity(:,j)=10.^((1./average).*sum(Block(n).spectra(:,j*average -average +1:j*average),2));
%          result.Block(n).attenuation(:,j)=Block(n).reference - (1./average).*sum(Block(n).spectra(:,j*average -average +1:j*average),2);
     end
   end
end
toc
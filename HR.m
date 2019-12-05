function[sig_out]=HR(data,skok,fvz) % vstup data-EKG, fvz-vzorkovací frekvence, skok- decimace odpovídající výpoètu spektrogramu (výstup prùbìh tepové frekvence)
delka=length(data);
%decim=1;
 decim=1/(100/fvz); %% decimace na fvz 100 Hz (pro urychlení)
 fvz=fvz/decim;
 data=resample(data,1,decim);
% filter=fir1(26,2/(fvz/2),'high');
% data=filtfilt(filter,1,data);
% data=abs(hilbert(data)).^2;
% skok1=200;
% prah(1:floor((length(data)-fvz*10)/skok1)+1)=single(0);
% 
% mm=median(single(data));
% s=std(single(data));
% %% stanovení adaptivního prahu pro QRS detekci
% 
% parfor(i=0:floor((length(data)-fvz*10-1)/skok1))
%     
%     m=max(data(1+i*skok1:1+skok1*i+fvz*10));
%     if(1/2*m>0.3*(mm+s))
%         
%     prah(i+1)=1/2*m;
%     else
%       prah(i+1)=0.3*(mm+s);  
%     end
% end
% 
% prah(prah==0)=prah(round((length(data)-fvz*10)/skok1));
% 
% 
% %% detekce QRS s adaptivním prahem
% 
% [vel,lok]=findpeaks(data);
% lok=single(lok);
% vel=single(vel);
% clear data;
% parfor(i=1:length(lok))
%     if(ceil(lok(i)/skok1)<length(prah))
%         if(vel(i)<prah(ceil(lok(i)/skok1)))
%             lok(i)=0;
%         end
%     else
%         if(vel(i)<prah(end))
%             lok(i)=0;
%         end
%     end
%     
%     
% end
% %clear data;
% %vel(lok==0)=[];
% lok(lok==0)=[];
% 
% 
% clear vel;
% clear prah;


lok=QRSdet_WAVE_5(fvz,data); %výpoèet pozic R kmitù pomocí detektoru s vlnkovou transformací
%% vytvoøení vektoru RR interval
RR1(1:length(lok)-1)=double(0);
i=int32(1);
lok=lok./fvz;
parfor(i=2:length(lok)) 
  
    RR1(i-1)= lok(i)-lok(i-1);
    
end

RR=double(RR1);

lok(1)=0;

RR=[RR RR(end)];
lok(end)=(delka/fvz-1/fvz);

%lok(RR<60/250)=60/250;
RR(RR<60/250)=60/250;

%RR=medfilt1(RR,9);
%% úprava výsledného signálu tepové frekvence
sig_out=resample(RR,double(lok),fvz/(skok/decim)); % pøevzorkování na shodnou fvz jako ostatní pøíznaky
sig_out=60./sig_out; % pøevod na TF  [/min]


% oøezání odlehlých hodnot
%sig_out(abs(sig_out)>350)=350;
sig_out(sig_out<0)=0;
% filtrace dolní propustí pro vyhlazení prùbìhu TF
%filter=fir1(1500,[1.5/(fvz/2)],'low');
%sig_out=filtfilt(filter,1,sig_out);
%sig_out=sig_out(1:round(length(data)/10));
end
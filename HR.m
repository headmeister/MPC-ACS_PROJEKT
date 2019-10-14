function[sig_out]=HR(data,skok,fvz) % vstup data-EKG, fvz-vzorkovac� frekvence, skok- decimace odpov�daj�c� v�po�tu spektrogramu (v�stup pr�b�h tepov� frekvence)

filter=fir1(1000,[8/(fvz/2) 20/(fvz/2)],'bandpass');
data=filtfilt(filter,1,data);
data=abs(hilbert(data)).^2;

m=median(data);
s=std(data);

prah=m+3*s;

%% detekce QRS
[vel,lok]=findpeaks(data); 

lok(vel<prah)=[];
lok=lok./fvz;

%% vytvo�en� vektoru RR interval
RR(1:length(lok)-1)=0;
for(i=2:length(lok)) 
RR(i-1)= lok(i)-lok(i-1);   

end
lok(1)=0;
RR=[RR(1) RR RR(end)];
lok=[lok (length(data)/fvz-1/fvz)];

%% �prava v�sledn�ho sign�lu tepov� frekvence
sig_out=resample(RR,lok,fvz/skok,'spline'); % p�evzorkov�n� na shodnou fvz jako ostatn� p��znaky
sig_out=60./sig_out; % p�evod na TF  [/min]


% o�ez�n� odlehl�ch hodnot
sig_out(abs(sig_out)>350)=350;
sig_out(sig_out<0)=0;
% filtrace doln� propust� pro vyhlazen� pr�b�hu TF
%filter=fir1(1500,[1.5/(fvz/2)],'low');
%sig_out=filtfilt(filter,1,sig_out);
%sig_out=sig_out(1:round(length(data)/10));
end
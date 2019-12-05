function[kval_v] = hodnoceniKvality(A,Acc)
load('sit.mat')
fis=readfis('FIS_std_fce.fis');
delka_sig=length(A);
skok=1000;
fvz=1000;

progress=0;
A=A-mean(A);
%% uprava dat akcelerometru (odstranìní SS složky + výpoèet celkové energie ve všech osách)

filter=fir1(500,[0.7/(100/2)],'high');

Acc(1,:)=filtfilt(filter,1,Acc(1,:));
Acc(2,:)=filtfilt(filter,1,Acc(2,:));
Acc(3,:)=filtfilt(filter,1,Acc(3,:));



vyk=sqrt(Acc(1,:).^2+Acc(2,:).^2+Acc(3,:).^2);
clear Acc;
%% výpoèet pøíznakù
accel_prum=skew(vyk,round(skok/(fvz/100)),1); %støední hodnota v èase
accel_std=skew(vyk,round(skok/(fvz/100)),2); % rozptyl v èase
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
clear Acc;

snr=SNR(A,skok,fvz); % odhad SNR
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

rate=HR(A,skok,fvz); % tepová frekvence
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);


power=spektrogram(A,skok); % vytvoøení spektrogramu
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
sikmost=skew(A,skok,3); % vytvoøení skew EKG
%clear A;
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
spicatost=skew(A,skok,4); % vytvoøení kurtosis EKG

drift=drift_rel(power,skok,fvz,0,3);% zjištìní relativní spektrální mocnosti driftu
progress=progress+5;
disp(['Progress: ' num2str(progress) '%']);
myo=drift_rel(power,skok,fvz,70,fvz/2); % zjištìní relativní spektrální mocnosti VF složek
progress=progress+5;
disp(['Progress: ' num2str(progress) '%']);
brum=drift_rel(power,skok,fvz,49,51); % zjištìní relativní spektrální mocnosti 50 Hz složky
progress=progress+5;
disp(['Progress: ' num2str(progress) '%']);




%% oøezání dat na stejnou délku
delka=min([length(accel_prum),length(accel_std),length(sikmost),length(spicatost),length(drift),length(myo),length(rate),length(brum)]); 


accel_prum=accel_prum(1:delka);
accel_std=accel_std(1:delka);
%accel_pow=(accel_pow-mean(accel_pow))/std(accel_pow);
sikmost=sikmost(1:delka);
%sikmost=(sikmost-mean(sikmost))/std(sikmost);
spicatost=spicatost(1:delka);
drift=drift(1:delka);
%drift=(drift-mean(drift))/std(drift);
myo=myo(1:delka);
%myo=(myo-mean(myo))/std(myo);
rate=rate(1:delka);
%rate=(rate-mean(rate))/std(rate);
brum=brum(1:delka);
%brum=(brum-mean(brum))/std(brum);
snr=snr(1:delka);
%snr=(snr-mean(snr))/std(snr);

% diference=diference(1:delka);
% diference=(diference-mean(diference))/std(diference);
clear delka;
disp('Vypocet priznaku hotov');


%% vybavovani

P(1:9,1:length(myo))=0;

P(1,:)=drift(1:end);
P(2,:)=myo(1:end);
P(3,:)=brum(1:end);
P(4,:)=rate(1:end);
P(5,:)=snr(1:end);
P(6,:)=sikmost(1:end);
P(7,:)=spicatost(1:end);
P(8,:)=accel_prum(1:end);
P(9,:)=accel_std(1:end);


P(isnan(P))=0; 

save('P3.mat','P');
progress=progress+5;
disp(['Progress: ' num2str(progress) '%']);

FF=sit(P,'useParallel','yes'); % vybavení sítì ze vstupu pøíznakù

progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

    kval=FF(1,:)+FF(2,:)*2+FF(3,:)*3; % vytvoøení signálu kvality
    %kval=[1 kval 1];
disp('vybaveno'); 
    
    
  
   matrix(1:length(kval)-2,1:3)=0;
  for(i=2:length(kval)-1)
      matrix(i-1,:)=kval(i-1:i+1);
      
      
  end  
  
kval(2:end-1)=evalfis(matrix,fis); % filtrace pomocí fis filtru
  
  
%kval=[pom(1)', kval', pom(end)'];
    
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
disp('FIS evaluace kompletní');    
 %kval= vec2ind(round(FF)); % návrat ke kvalitì
 kval_v=(resample(kval,round(delka_sig/length(FF(1,:))),1)); % pøevzorkování zpìt

kval_v(kval_v>3)=3;
kval_v(kval_v<1)=1;


kval_v(kval_v==0)=1;
kval_v=round(kval_v);

if(length(kval_v)>delka_sig) % úprava délky
    ppp=length(kval_v)-delka_sig;
    kval_v=kval_v(1:end-ppp);
elseif (length(kval_v)<delka_sig)
    ppp=delka_sig-length(kval_v);
    ppp=kval_v(end)*ones(1,ppp);
    kval_v=[kval_v ppp];
end
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
disp('Hotovo');   
end
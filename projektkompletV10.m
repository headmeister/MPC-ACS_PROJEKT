clc
clear
close all
fis=readfis('FIS_std_fce.fis');
[head,record]=edfread('100001.EDF');
load('100001_ann_quality_train.mat');
an=anotace;
load('105001_ann_quality_train.mat');
an=[an anotace];
load('111001_ann_quality_train.mat');
anotace=int8([an anotace]);
%plot(record(1,:));
A1=(record(1,1:68120000));
Acc=record(2:4,1:6812000);
clear record;

[head,record]=edfread('105001.EDF');
A1=[A1 record(1,20000000:80000000)];
Acc=[Acc record(2:4,2000000:8000000)];
clear record;
[head,record]=edfread('111001.EDF');
A=[A1 record(1,1:30000000)];
Acc=[Acc record(2:4,1:3000000)];
clear A1;
clear an;



skok=500;
fvz=1000;
pocet_vzorku=100000; % poèet vzorkù od jednotlivých kategorií v množinì pro uèení
progress=0;
%% výpoèet pøíznakù
accel_pow=accel(Acc,fvz/skok,100); % odhad úrovnì akcelerace akcelerometru (bez SS složky)
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
sikmost=skew(A,skok); % vytvoøení skew
%clear A;
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
%spicatost=kurt(A,skok); % vytvoøení skew

drift=drift_rel(power,skok,fvz,0,2);% zjištìní relativní spektrální mocnosti driftu
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
myo=drift_rel(power,skok,fvz,100,fvz/2); % zjištìní relativní spektrální mocnosti VF složek
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
brum=drift_rel(power,skok,fvz,49,51); % zjištìní relativní spektrální mocnosti 50 Hz složky
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

% diference=abs(diff(rate));

delka=min([length(accel_pow),length(sikmost),length(drift),length(myo),length(rate),length(brum),length(diference)]);

accel_pow=accel_pow(1:delka);
accel_pow=(accel_pow-mean(accel_pow))/std(accel_pow);
sikmost=sikmost(1:delka);
sikmost=(sikmost-mean(sikmost))/std(sikmost);
%spicatost=spicatost(1:delka);
drift=drift(1:delka);
drift=(drift-mean(drift))/std(drift);
myo=myo(1:delka);
myo=(myo-mean(myo))/std(myo);
rate=rate(1:delka);
rate=(rate-mean(rate))/std(rate);
brum=brum(1:delka);
brum=(brum-mean(brum))/std(brum);
snr=snr(1:delka);
snr=(snr-mean(snr))/std(snr);

% diference=diference(1:delka);
% diference=(diference-mean(diference))/std(diference);
disp('Vypocet priznaku hotov');

%% neuro

a_i=1;
b_i=1;
c_i=1;
ind=1;
anotace_pom=round(resample(double(anotace),1,skok));
anotace_pom=anotace_pom(1:length(myo));
pozice1=int64(find(anotace_pom==1));
pozice2=int64(find(anotace_pom==2));
pozice3=int64(find(anotace_pom==3));

maxima(1)=length(pozice1);
maxima(2)=length(pozice2);
maxima(3)=length(pozice3);

priznak=zeros(7,pocet_vzorku*3);


an=logical(zeros(3,pocet_vzorku*3));






for(i=int64(1:3))
    
    if(i==1)
        if(pocet_vzorku>=maxima(i))
        Q=maxima(i);
        else
         Q=pocet_vzorku; 
         
        end
        index1=zeros(2,Q);
    index1(2,:)=randperm(length(pozice1),Q); 
    index1(2,:)=pozice1(index1(2,1:Q));
    index1(1,:)=i;
    elseif(i==2)
        if(pocet_vzorku>=maxima(i))
        Q=maxima(i);
        else
         Q=pocet_vzorku;   
        end
        index2=zeros(2,Q);
        index2(2,:)=randperm(length(pozice2),Q);    
        index2(2,:)=pozice2(index2(2,1:Q));
        index2(1,:)=i;
    else
        if(pocet_vzorku>=maxima(i))
        Q=maxima(i);
        else
         Q=pocet_vzorku;   
        end
            index3=zeros(2,Q);
          index3(2,:)=randperm(length(pozice3),Q);    
          index3(2,:)=pozice3(index3(2,1:Q));
          index3(1,:)=i;
    end
end

index=[index1 index2 index3];

for(i=int64(1:length(index(1,:))))
        priznak(1,i)=drift(index(2,i));
        priznak(2,i)=myo(index(2,i));
        priznak(3,i)=accel_pow(index(2,i));
        priznak(4,i)=brum(index(2,i));
        priznak(5,i)=rate(index(2,i));
        priznak(6,i)=snr(index(2,i));
        priznak(7,i)=sikmost(index(2,i));
       % priznak(8,i)=diference(index(2,i));
        an(index(1,i),i)=1;
end

priznak(:,length(index)+1:end)=[];
an(:,length(index)+1:end)=[];

priznak(isnan(priznak))=0;
%% uceni/vybavovani
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

sit=patternnet(8); % vytvoøení sítì s n skrytými neurony

sit=train(sit,priznak,an,'useParallel','yes'); % uèení sítì na trénovací množinì
FF=sit([drift;myo;accel_pow;brum;rate;snr;sikmost],'useParallel','yes'); % vybavení sítì ze vstupu pøíznakù

progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

    pom=FF(1,:)+FF(2,:)*2+FF(3,:)*3;
disp('vybaveno'); 
    
    
  matrix(1:length(pom)-2,1:3)=0;
  for(i=2:length(pom)-1)
      matrix(i-1,:)=pom(i-1:i+1);
  end  
   
    
   kval=evalfis(matrix,fis);
   kval=[pom(1)', kval', pom(end)'];
    
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
disp('FIS evaluace kompletní');    
 %kval= vec2ind(round(FF)); % návrat ke kvalitì
 kval=round(resample(kval,round(length(anotace)/length(FF(1,:))),1)); % pøevzorkování zpìt
% 
%kval=medfilt1(round(kval),3*skok); % spojení do delších úsekù
kval(kval==0)=1;
kval=round(kval);

if(length(kval)>length(anotace))
    ppp=length(kval)-length(anotace);
    kval=kval(1:end-ppp);
elseif (length(kval)<length(anotace))
    ppp=length(anotace)-length(kval);
    ppp=kval(end)*ones(1,ppp);
    kval=[kval ppp];
    
end


 [uspesnost, uspescnost_cn]=hodnoceni(kval,anotace);

[SE, PP, F1, cov]=F1score(anotace,kval);
disp(['SE: ' num2str(SE,3) ' %' ]);
disp(['PP: ' num2str(PP,3) ' %']);
disp(['F1: ' num2str(F1,3) ', coverage: ' num2str(cov,3) ' %']);
% %% analýza korelace
% Korelace=[];
% anotace_short= vec2ind(an);
% for i=1:length(priznak(:,1))
% Korelace(i,1)=corr(anotace_short',priznak(i,:)','Type','Pearson');
% Korelace(i,2)=corr(anotace_short',priznak(i,:)','Type','Spearman');
% %Korelace(i,3)=corr(anotace_short',priznak(i,:)','Type','Kendall');
% end


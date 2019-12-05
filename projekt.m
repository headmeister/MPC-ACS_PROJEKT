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
pocet_vzorku=200000; % po�et vzork� od jednotliv�ch kategori� v mno�in� pro u�en�
progress=0;
%% v�po�et p��znak�
accel_pow=accel(Acc,fvz/skok,100); % odhad �rovn� akcelerace akcelerometru (bez SS slo�ky)
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
clear Acc;

snr=SNR(A,skok,fvz); % odhad SNR
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

rate=HR(A,skok,fvz); % tepov� frekvence
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);


power=spektrogram(A,skok); % vytvo�en� spektrogramu
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
sikmost=skew(A,skok); % vytvo�en� skew
%clear A;
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
%spicatost=kurt(A,skok); % vytvo�en� skew

drift=drift_rel(power,skok,fvz,0,2);% zji�t�n� relativn� spektr�ln� mocnosti driftu
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
myo=drift_rel(power,skok,fvz,100,fvz/2); % zji�t�n� relativn� spektr�ln� mocnosti VF slo�ek
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
brum=drift_rel(power,skok,fvz,49,51); % zji�t�n� relativn� spektr�ln� mocnosti 50 Hz slo�ky
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);


delka=min([length(accel_pow),length(sikmost),length(drift),length(myo),length(rate),length(brum)]);

accel_pow=accel_pow(1:delka);
sikmost=sikmost(1:delka);
%spicatost=spicatost(1:delka);
drift=drift(1:delka);
myo=myo(1:delka);
rate=rate(1:delka);
brum=brum(1:delka);
snr=snr(1:delka);

disp('Vypocet priznaku hotov');

% P(1:6,1:length(myo))=0;
% 
% P(1,:)=accel_pow;
% P(2,:)=drift;
% P(3,:)=myo;
% P(4,:)=rate;
% P(5,:)=brum;
% %P(6,:)=snr;




% subplot(311);
% plot(A);
% subplot(312);
% plot(drift);
% subplot(313);
% plot(myo);





%% jednoduch� hodnocen� kvality na z�klad� prah� a dvou p��znak� a zhodnocen� �sp�nosti
%kvalita(1:length(power))=0;
% for(i=1:length(power))
%
%
%     if((drift(i)>0.6&&drift(i)<0.95) || (myo(i)>0.015&& myo(i)<0.02))
%         kvalita(i)=2;
%     elseif(drift(i)>=0.95 ||myo(i)>=0.02)
%         kvalita(i)=3;
%     else
%         kvalita(i)=1;
%     end
%     a=kvalita(i);
%     aa=anotace(i);
%     if(kvalita(i)==1)
%         uspesnost(1,1)=uspesnost(1,1)+(kvalita(i)==anotace(i));
%
%     elseif(kvalita(i)==2)
%         uspesnost(2,1)=uspesnost(2,1)+(kvalita(i)==anotace(i));
%     else
%         uspesnost(3,1)=uspesnost(3,1)+(kvalita(i)==anotace(i));
%     end
% end







%% neuro

%an=anotace(1:10000);
%
% an(1,:)=anotace==1;
% an(2,:)=anotace==2;
% an(3,:)=anotace==3;


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

%priznak(1:6,1:pocet_vzorku*3)=0;
priznak=zeros(6,pocet_vzorku*3);


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
        %priznak(7,i)=spicatost(index(2,i));
        an(index(1,i),i)=1;
end

priznak(:,length(index)+1:end)=[];
an(:,length(index)+1:end)=[];
% while((a_i<=pocet_vzorku || b_i<=pocet_vzorku || c_i<=pocet_vzorku)&&(a_i<sum(anotace==1)&&b_i<sum(anotace==2)&&c_i<sum(anotace==3))) % v�b�r u��c� mno�iny dat
%     if(anotace(index)==1 && a_i<=pocet_vzorku)
%         priznak(1,ind)=drift(index);
%         priznak(2,ind)=myo(index);
%         priznak(3,ind)=accel_pow(index);
%         priznak(4,ind)=brum(index);
%         priznak(5,ind)=rate(index);
%         priznak(6,ind)=snr(index);
%         an(1,ind)=1;
%         ind=ind+1;
%         a_i=a_i+1;
%     end
%     
%     if(anotace(index)==2 && b_i<=pocet_vzorku)
%         priznak(1,ind)=drift(index);
%         priznak(2,ind)=myo(index);
%         priznak(3,ind)=accel_pow(index);
%         priznak(4,ind)=brum(index);
%         priznak(5,ind)=rate(index);
%         priznak(6,ind)=snr(index);
%         an(2,ind)=1;
%         ind=ind+1;
%         b_i=b_i+1;
%     end
%     
%     if(anotace(index)==3 && c_i<=pocet_vzorku)
%         priznak(1,ind)=drift(index);
%         priznak(2,ind)=myo(index);
%         priznak(3,ind)=accel_pow(index);
%         priznak(4,ind)=brum(index);
%         priznak(5,ind)=rate(index);
%         priznak(6,ind)=snr(index);
%         an(3,ind)=1;
%         ind=ind+1;
%         c_i=c_i+1;
%     end
%     
%     
%            
% end
priznak(isnan(priznak))=0;
%% uceni/vybavovani
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

sit=patternnet(5); % vytvo�en� s�t� s 7 skryt�mi vrstvami

sit=train(sit,priznak,an,'useParallel','yes'); % u�en� s�t� na tr�novac� mno�in�
FF=sit([drift;myo;accel_pow;brum;rate;snr;sikmost],'useParallel','yes'); % vybaven� s�t� ze vstupu p��znak�

progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);

    pom=FF(1,:)+FF(2,:)*2+FF(3,:)*3;
disp('vybaveno'); 
   
    %h = waitbar(0, 'Please wait ...');
    
    
  matrix(1:length(pom)-2,1:3)=0;
  for(i=2:length(pom)-1)
      matrix(i-1,:)=pom(i-1:i+1);
  end  
   
    
   kval=evalfis(matrix,fis);
   kval=[pom(1)', kval', pom(end)'];
    
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
disp('FIS evaluace kompletn�');    
 %kval= vec2ind(round(FF)); % n�vrat ke kvalit�
 kval=round(resample(kval,round(length(anotace)/length(FF(1,:))),1)); % p�evzorkov�n� zp�t
% 
%kval=medfilt1(round(kval),3*skok); % spojen� do del��ch �sek�
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

%% vykreslen� v�sledn�ho hodnocen� kvality a anotace
% subplot(211);
% plot(anotace);
% subplot(212);
% plot(kval);

%% anal�za korelace
Korelace=[];
anotace_short= vec2ind(an);
for i=1:length(priznak(:,1))
Korelace(i,1)=corr(anotace_short',priznak(i,:)','Type','Pearson');
Korelace(i,2)=corr(anotace_short',priznak(i,:)','Type','Spearman');
%Korelace(i,3)=corr(anotace_short',priznak(i,:)','Type','Kendall');
end

 function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end

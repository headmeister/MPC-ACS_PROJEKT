clc
clear
close all
[head,record]=edfread('100001.EDF');
load('100001_ann_quality_train.mat');
%plot(record(1,:));
A=(record(1,1:68120000));
skok=500;
fvz=1000;
pocet_vzorku=100000; % po�et vzork� od jednotliv�ch kategori� v mno�in� pro u�en�
%% v�po�et p��znak�
accel_pow=accel(record(2:4,1:6812000),fvz/skok,100); % v�kon akcelerometru (bez SS slo�ky)
power=spektrogram(A,skok); % vytvo�en� spektrogramu
sikmost=skew(A,skok); % vytvo�en� skew
spicatost=kurt(A,skok); % vytvo�en� skew

drift=drift_rel(power,skok,fvz,0,2);% zji�t�n� relativn� spektr�ln� mocnosti driftu
myo=drift_rel(power,skok,fvz,100,fvz/2); % zji�t�n� relativn� spektr�ln� mocnosti VF slo�ek
rate=HR(A,skok,fvz); % tepov� frekvence
brum=drift_rel(power,skok,fvz,49,51); % zji�t�n� relativn� spektr�ln� mocnosti 50 Hz slo�ky
snr=SNR(A,skok,fvz); % odhad SNR




P(1:6,1:length(myo))=0;

P(1,:)=accel_pow;
P(2,:)=drift;
P(3,:)=myo;
P(4,:)=rate;
P(5,:)=brum;
P(6,:)=snr;




% subplot(311);
% plot(A);
% subplot(312);
% plot(drift);
% subplot(313);
% plot(myo);



uspesnost(1:3,1)=0;

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
anotace_pom=round(resample(anotace,1,skok));

pozice1=int64(find(anotace_pom==1));
pozice2=int64(find(anotace_pom==2));
pozice3=int64(find(anotace_pom==3));

maxima(1)=length(pozice1);
maxima(2)=length(pozice2);
maxima(3)=length(pozice3);

%priznak(1:6,1:pocet_vzorku*3)=0;
priznak=zeros(6,pocet_vzorku*3);


an=logical(zeros(3,pocet_vzorku*3));






for(i=1:3)
    
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


for(i=1:length(index(1,:)))
        priznak(1,i)=drift(index(2,i));
        priznak(2,i)=myo(index(2,i));
        priznak(3,i)=accel_pow(index(2,i));
        priznak(4,i)=brum(index(2,i));
        priznak(5,i)=rate(index(2,i));
        priznak(6,i)=snr(index(2,i));
        priznak(7,i)=sikmost(index(2,i));
        priznak(8,i)=spicatost(index(2,i));
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


%% uceni/vybavovani
sit=patternnet(50); % vytvo�en� s�t� s 25 skryt�mi vrstvami

sit=train(sit,priznak,an,'useParallel','no'); % u�en� s�t� na tr�novac� mno�in�
FF=sit([drift; myo;accel_pow;brum;rate;snr;sikmost;spicatost],'useParallel','no'); % vybaven� s�t� ze vstupu p��znak�



kval= vec2ind(round(FF)); % n�vrat ke kvalit�
kval=round(resample(kval,skok,1)); % p�evzorkov�n� zp�t
kval=medfilt1(round(kval),3*skok); % spojen� do del��ch �sek�


%% ur�en� �sp�nosti
for(i=1:length(kval))
    if(kval(i)==1)
        uspesnost(1,1)=uspesnost(1,1)+(kval(i)==anotace(i));
        
    elseif(kval(i)==2)
        uspesnost(2,1)=uspesnost(2,1)+(kval(i)==anotace(i));
    else
        uspesnost(3,1)=uspesnost(3,1)+(kval(i)==anotace(i));
    end
    
end


uspesnost(1)=uspesnost(1)/sum(anotace==1)*100;

uspesnost(2)=uspesnost(2)/sum(anotace==2)*100;

uspesnost(3)=uspesnost(3)/sum(anotace==3)*100;

uspesnost_cn= sum(kval==anotace)/length(anotace)*100;
%% vykreslen� v�sledn�ho hodnocen� kvality a anotace
subplot(211);
plot(anotace);
subplot(212);
plot(kval);
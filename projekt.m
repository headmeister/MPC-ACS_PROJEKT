clc
clear 
close all
[head,record]=edfread('100001.EDF');
load('100001_ann_quality_train.mat');
%plot(record(1,:));
A=record(1,1:68120000);
skok=1000;
fvz=1000;
%% v�po�et p��znak�
accel_pow=accel(record(2:4,1:6812000),1000,100); % v�kon akcelerometru (bez SS slo�ky)
power=spektrogram(A,skok); % vytvo�en� spektrogramu
drift=drift_rel(power,skok,fvz,0,2);% zji�t�n� relativn� spektr�ln� mocnosti driftu
myo=drift_rel(power,skok,fvz,100,fvz/2); % zji�t�n� relativn� spektr�ln� mocnosti VF slo�ek
rate=HR(A,1000); % tepov� frekvence
brum=drift_rel(power,skok,fvz,49,51); % zji�t�n� relativn� spektr�ln� mocnosti 50 Hz slo�ky

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
pocet_vzorku=100000; % po�et vzork� od jednotliv�ch kategori� v mno�in� pro u�en�

a_i=1;
b_i=1;
c_i=1;
index=1;
ind=1;

priznak(1:3,1:pocet_vzorku*3)=0;
an(1:3,1:pocet_vzorku*3)=0;






while((a_i<=pocet_vzorku || b_i<=pocet_vzorku || c_i<=pocet_vzorku)&&index<length(drift)) % v�b�r u��c� mno�iny dat
if(anotace(index)==1 && a_i<=pocet_vzorku)
priznak(1,ind)=drift(index);
priznak(2,ind)=myo(index);
priznak(3,ind)=accel_pow(index);
priznak(4,ind)=brum(index);
priznak(5,ind)=rate(index);
an(1,ind)=1;
ind=ind+1;
a_i=a_i+1;
end

if(anotace(index)==2 && b_i<=pocet_vzorku)
priznak(1,ind)=drift(index);
priznak(2,ind)=myo(index);
priznak(3,ind)=accel_pow(index);
priznak(4,ind)=brum(index);
priznak(5,ind)=rate(index);
an(2,ind)=1;
ind=ind+1;
b_i=b_i+1;
end

if(anotace(index)==3 && c_i<=pocet_vzorku)
priznak(1,ind)=drift(index);
priznak(2,ind)=myo(index);
priznak(3,ind)=accel_pow(index);
priznak(4,ind)=brum(index);
priznak(5,ind)=rate(index);
an(3,ind)=1;
ind=ind+1;
c_i=c_i+1;
end
index=index+1;
end



sit=patternnet(10); % vytvo�en� s�t� s 10 skryt�mi vrstvami

sit=train(sit,priznak,an); % u�en� s�t� na tr�novac� mno�in�
FF=sit([drift; myo;accel_pow;brum;rate]); % vybaven� s�t� ze vstupu p��znak�



kval= vec2ind(round(FF)); % n�vrat ke kvalit�
kval=medfilt1(kval,3000); % spojen� do del��ch �sek�

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
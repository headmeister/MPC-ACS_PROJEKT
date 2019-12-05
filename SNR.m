function[SNR]=SNR(data,skok1,fvz) %Funkce pro odhad SNR na z�klad� filtrace sign�lu pomoc� SWT
% data- sign�l EKG
% vstupy: skok1= okno resp decimace vstupn�ho sign�lu
% fvz- vzorkovac� fce sign�lu
decim=round(1/(500/fvz));
fvz=fvz/decim;
data=resample(data,1,decim);
filter=fir1(500,[2/(fvz/2)],'high');
data=filtfilt(filter,1,data);

vlnka = 'sym4';
stup_rozk = 4;
fragmentace=8;
for(i=0:2^stup_rozk*fragmentace)
    if(mod(length(data)+i,2^stup_rozk*fragmentace)==0)
        
        pom=i;
        break
    end
end

data(end:end+pom)=0;


fragment=length(data)/fragmentace;
SNR=data;
for(i=1:fragmentace)
    SWC = single(swt(data(1+(i-1)*fragment:i*fragment), stup_rozk, vlnka));
    prahy=single(SWC);
    %% prahovani
    skok=500000;
    konstanta=0.0890;
    for(j=1:skok:length(SWC(1,:))-skok)
        prahy(1:stup_rozk,j:j+skok)=repmat(median(abs(SWC(1:stup_rozk,j:j+skok)),2)/konstanta,1,skok+1);
    end
    prahy(1:stup_rozk,j:end)=repmat(prahy(1:stup_rozk,j),1,length(prahy)-j+1);
    
    SWC(1:stup_rozk,:)=(SWC(1:stup_rozk,:)-prahy(1:stup_rozk,:).^2./SWC(1:stup_rozk,:)).*(abs(SWC(1:stup_rozk,:))>prahy(1:stup_rozk,:));
    
    clear prahy;
    
    SWC=double(SWC);
    
    data2 = iswt(SWC, vlnka);
    clear SWC;
    
    
    noise=data(1+(i-1)*fragment:i*fragment)-data2;
    SNR(1+(i-1)*fragment:i*fragment)=10*log10(((data2).^2./(noise).^2));
    
end
SNR(SNR>10000)=10000;
SNR(end-pom+1:end)=[];
SNR=resample(SNR,1,skok1/decim);
% filter=fir1(1500,[2/(fvz/2)],'low');
% SNR=filtfilt(filter,1,SNR);


end
function[SNR]=SNR(data,skok,fvz) %Funkce pro odhad SNR na základì odhadu spektrálních složek užiteèného signálu a šumu( 50Hz, 0-2Hz,70-fvz Hz)

filter=fir1(1500,[2/(fvz/2) 70/(fvz/2)],'bandpass');
data2=filtfilt(filter,1,data);

filter=fir1(1500,[49/(fvz/2) 51/(fvz/2)],'stop');
data2=filtfilt(filter,1,data2);


noise=data-data2;

SNR=10*log((data2.^2)./(noise.^2));
SNR=resample(SNR,1,skok);
% filter=fir1(1500,[2/(fvz/2)],'low');
% SNR=filtfilt(filter,1,SNR);


end
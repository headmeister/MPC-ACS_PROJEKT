function [vystup]=TKEO(signal,skok)
% signal == vstupní signál
% skok   == skok pro pøevzrokování
% vystup == TKEO
vystup=zeros(1,length(signal));
for i=1:length(signal)
    if i==1
        vystup(i)=signal(i)^2-signal(i+1)^2;
    elseif i==length(signal)
        vystup(i)=signal(i)^2-signal(i-1)^2;
    else
        vystup(i)=signal(i)^2-signal(i+1)*signal(i-1);
    end
end
vystup=resample(vystup,1,skok);
end
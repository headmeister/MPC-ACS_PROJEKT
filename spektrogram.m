function[power]=spektrogram(data,skok) %vytvaøí spektrogram bez pøekrývání z vektoru data po úsecích o délce skok
n=1;
power(1:length(data))=0;
while((n+skok)<length(data))
    power(n:n+skok)=abs(fft(data(n:n+skok))).^2;
    power(n)=0.001; % odstranìní SS složky má pramalý vliv na hodnocení kvality EKG
    n=n+skok;
end
power(n:end)=abs(fft(data(n:end))).^2;

end
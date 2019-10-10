function[power]=spektrogram(data,skok) %vytva�� spektrogram bez p�ekr�v�n� z vektoru data po �sec�ch o d�lce skok
n=1;
power(1:length(data))=0;
while((n+skok)<length(data))
    power(n:n+skok)=abs(fft(data(n:n+skok))).^2;
    power(n)=0.001; % odstran�n� SS slo�ky m� pramal� vliv na hodnocen� kvality EKG
    n=n+skok;
end
power(n:end)=abs(fft(data(n:end))).^2;

end
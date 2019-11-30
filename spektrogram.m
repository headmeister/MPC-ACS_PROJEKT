function[power]=spektrogram(data,skok) %vytva�� spektrogram bez p�ekr�v�n� z vektoru data po �sec�ch o d�lce skok
% n=1;
% power(1:length(data))=0;
% while((n+skok)<length(data))
%     pom=data(n:n+skok-1);
%     pom=pom.*gausswin(skok)';
%     power(n:n+skok-1)=abs(fft(pom)).^2;
%     power(n)=0.001; % odstran�n� SS slo�ky m� pramal� vliv na hodnocen� kvality EKG
%     n=n+skok;
% end
% if(n<length(data))
% pom=data(n:end);
% pom=pom.*gausswin(length(data)-n+1)';
% power(n:end)=abs(fft(pom)).^2;
% end
power=abs(spectrogram(data,skok,0,skok)).^2;

end
function[power]=spektrogram(data,skok) %vytvaøí spektrogram bez pøekrývání z vektoru data po úsecích o délce skok
% n=1;
% power(1:length(data))=0;
% while((n+skok)<length(data))
%     pom=data(n:n+skok-1);
%     pom=pom.*gausswin(skok)';
%     power(n:n+skok-1)=abs(fft(pom)).^2;
%     power(n)=0.001; % odstranìní SS složky má pramalý vliv na hodnocení kvality EKG
%     n=n+skok;
% end
% if(n<length(data))
% pom=data(n:end);
% pom=pom.*gausswin(length(data)-n+1)';
% power(n:end)=abs(fft(pom)).^2;
% end
power=abs(spectrogram(data,skok,0,skok)).^2;

end
function[skew]=skew(data,skok) %vytva�� skew bez p�ekr�v�n� z vektoru data po �sec�ch o d�lce skok
n=1;

ind=1;
for(i=0:skok)
    if(mod(length(data)+i,skok)==0)
        
        pom=i;
        break
    end
end

data(end:end+pom)=0;


data=reshape(data,skok,length(data)/skok);
skew(1:length(data(1,:)))=0;

skew(1:end-1)=skewness(data(:,1:end-1));
pom=data(end,:);
skew(end)=skewness(pom(pom~=0));


% while((n+skok)<length(data))
%     pom=data(n:n+skok-1);
%     skew(ind)=skewness(pom);
%     ind=ind+1;
%     n=n+skok;
% end
% 
% if(n<length(data))
% pom=data(n:end);
% skew(ind)=skewness(pom);
% end
end
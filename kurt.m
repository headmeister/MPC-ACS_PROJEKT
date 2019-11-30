function[skew1]=kurt(data,skok) %vytvaøí skew bez pøekrývání z vektoru data po úsecích o délce skok
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
skew1(1:length(data(1,:)))=0;

skew1(1:end-1)=kurtosis(data(:,1:end-1));
pom=data(end,:);
skew1(end)=kurtosis(pom(pom~=0));



end
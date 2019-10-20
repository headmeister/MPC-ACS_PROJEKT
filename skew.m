function[skew]=skew(data,skok) %vytvaøí skew bez pøekrývání z vektoru data po úsecích o délce skok
n=1;
skew(1:round(length(data)/skok))=0;
ind=1;
while((n+skok)<length(data))
    pom=data(n:n+skok-1);
    skew(ind)=skewness(pom);
    ind=ind+1;
    n=n+skok;
end

if(n<length(data))
pom=data(n:end);
skew(ind)=skewness(pom);
end
end
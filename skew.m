function[skew]=skew(data,skok) %vytva�� skew bez p�ekr�v�n� z vektoru data po �sec�ch o d�lce skok
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
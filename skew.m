function[skew1]=skew(data,skok,param) % v�stup statistick� momenty v okn� bez p�ekr�v�n� z vektoru data po �sec�ch o d�lce skok
%vstup- data (sign�l)
%skok- delka okna (stupe� decimace)
%param - statistick� moment (1-pr�m�r,2-std, 3-skew, 4- kurtosis)
for(i=0:skok)
    if(mod(length(data)+i,skok)==0)
        
        pom=i;
        break
    end
end

data(end:end+pom)=0;


data=reshape(data,skok,length(data)/skok);
skew1(1:length(data(1,:)))=0;
if(param==3)
skew1(1:end-1)=skewness(data(:,1:end-1));
pom=data(end,:);
skew1(end)=skewness(pom(pom~=0));
elseif (param==4)
skew1(1:end-1)=kurtosis(data(:,1:end-1));
pom=data(end,:);
skew1(end)=kurtosis(pom(pom~=0));  

elseif (param==1)
skew1(1:end-1)=mean(data(:,1:end-1));
pom=data(end,:);
skew1(end)=mean(pom(pom~=0)); 

elseif (param==2)
skew1(1:end-1)=std(data(:,1:end-1));
pom=data(end,:);
skew1(end)=std(pom(pom~=0)); 
    
end


end
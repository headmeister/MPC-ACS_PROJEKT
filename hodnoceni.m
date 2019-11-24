function[uspesnost,uspesnost_cn]=hodnoceni(kval,anotace)% v�po�et procen-
%tu�ln� shody dvou vektor� kval- vektor z�skan� hodnocen�m sign�lu
%                         anotace - vektor obdr�en� jako hodnocen� expertem
%                         uspesnost - shoda v jednotlivych kategori�ch
%                         uspesnost_cn - Hodnocen� celkov� shody
uspesnost(1:3,1)=0;
for(i=1:length(kval))
    if(kval(i)==1)
        uspesnost(1,1)=uspesnost(1,1)+(kval(i)==anotace(i));
        
    elseif(kval(i)==2)
        uspesnost(2,1)=uspesnost(2,1)+(kval(i)==anotace(i));
    else
        uspesnost(3,1)=uspesnost(3,1)+(kval(i)==anotace(i));
    end
    
end


uspesnost(1)=uspesnost(1)/sum(anotace==1)*100;

uspesnost(2)=uspesnost(2)/sum(anotace==2)*100;

uspesnost(3)=uspesnost(3)/sum(anotace==3)*100;

uspesnost_cn= sum(kval==anotace)/length(anotace)*100;
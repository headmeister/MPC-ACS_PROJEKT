function[uspesnost,uspesnost_cn]=hodnoceni(kval,anotace)% výpoèet procen-
%tuální shody dvou vektorù kval- vektor získaný hodnocením signálu
%                         anotace - vektor obdržený jako hodnocení expertem
%                         uspesnost - shoda v jednotlivych kategoriích
%                         uspesnost_cn - Hodnocení celkové shody
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
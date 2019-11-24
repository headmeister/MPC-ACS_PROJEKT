% vstupy:   anotace - vektor anotac� kvality EKG 
%           stanovenaKvalita - vektor algoritmem stanoven� kvality EKG pro ka�d� vzorek
% v�stupy:  SE - senzitivita
%           PP - pozitivn� prediktivita
%           F1 - pr�m�rn� F1 sk�re (pr�m�r F1 v�ech 3 t��d kvality)
%           cov - coverage, procentu�ln� vyj�d�en� pom�r spr�vn� za�azen�ch vzork�

function [SE, PP, F1, cov]=F1score(anotace, stanovenaKvalita)
    
    % kontrola shodn� velikosti vstupn�ch vektor�
    if length(anotace)~=length(stanovenaKvalita)
        f = warndlg('Vektory "anotace" a "stanovenaKvalita" nemaj� stejnou velikost.');
        return
    end
    
    L=length(anotace);

    % stanoven� po�tu vzork� v ka�d� skupin� kvality
    [jedna, dva, tri]=jdt(L,anotace);
    celkem=jedna+dva+tri;
    anotaceZastoupeni=[jedna, dva, tri];            % po�et vzork� v ka�d� skupin� kvality
    clear pom jedna dva tri

    % ov��en�, �e vektor "anotace" neobsahuje nedovolen� hodnoty
    if celkem~=length(anotace)
        f = warndlg('Vektor "anotace" obsahuje nedovolen� hodnoty (povolena pouze ��sla 1, 2, 3).');
        return
    end
    
    % ov��en�, �e vektor "stanovenaKvalita" neobsahuje nedovolen� hodnoty
    if sum(stanovenaKvalita>3)>0 || sum(stanovenaKvalita<1)>0 || isreal(stanovenaKvalita)==0
        f = warndlg('Vektor "stanovenaKvalita" obsahuje nedovolen� hodnoty (povolena pouze ��sla 1, 2, 3).');
        return
    end
    
    podobnost=stanovenaKvalita;         
    podobnost(anotace~=stanovenaKvalita)=4;              % na pozice, kde se stanoven� kvalita neshoduje s anotac� se d� 4
    
    [jedna, dva, tri]=jdt(L,podobnost);
    TP=[jedna, dva, tri];            % po�et vzork� v ka�d� skupin� kvality
    clear pom jedna dva tri
    
    SE=TP./anotaceZastoupeni*100;% procentu�ln� zastoupen� spr�vn� za�azen�ch vzork� v ka�d� skupin�
    
    [jedna, dva, tri]=jdt(L,stanovenaKvalita);
    kvalitaZastoupeni=[jedna, dva, tri];
    
    % v�po�et F1 sk�re
    FP=kvalitaZastoupeni-TP;     % po�et fale�n� pozitivn�ch vzork� v ka�d� skupin�
    PP=TP./(TP+FP)*100; % pozitivn� prediktivita
    F1dilci=1./((1./SE+1./PP)./2);
    F1=nanmean(F1dilci);
    
    % coverage
    cov=((sum(TP))/length(anotace))*100;
end

function [jedna, dva, tri]=jdt(L,v)
    pom=zeros(1,L);
    pom(v==1)=1;
    jedna=sum(pom);
    pom=zeros(1,L);
    pom(v==2)=1;
    dva=sum(pom);
    pom=zeros(1,L);
    pom(v==3)=1;
    tri=sum(pom);
end
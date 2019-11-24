% vstupy:   anotace - vektor anotací kvality EKG 
%           stanovenaKvalita - vektor algoritmem stanovené kvality EKG pro každý vzorek
% výstupy:  SE - senzitivita
%           PP - pozitivní prediktivita
%           F1 - prùmìrné F1 skóre (prùmìr F1 všech 3 tøíd kvality)
%           cov - coverage, procentuálnì vyjádøený pomìr správnì zaøazených vzorkù

function [SE, PP, F1, cov]=F1score(anotace, stanovenaKvalita)
    
    % kontrola shodné velikosti vstupních vektorù
    if length(anotace)~=length(stanovenaKvalita)
        f = warndlg('Vektory "anotace" a "stanovenaKvalita" nemají stejnou velikost.');
        return
    end
    
    L=length(anotace);

    % stanovení poètu vzorkù v každé skupinì kvality
    [jedna, dva, tri]=jdt(L,anotace);
    celkem=jedna+dva+tri;
    anotaceZastoupeni=[jedna, dva, tri];            % poèet vzorkù v každé skupinì kvality
    clear pom jedna dva tri

    % ovìøení, že vektor "anotace" neobsahuje nedovolené hodnoty
    if celkem~=length(anotace)
        f = warndlg('Vektor "anotace" obsahuje nedovolené hodnoty (povolena pouze èísla 1, 2, 3).');
        return
    end
    
    % ovìøení, že vektor "stanovenaKvalita" neobsahuje nedovolené hodnoty
    if sum(stanovenaKvalita>3)>0 || sum(stanovenaKvalita<1)>0 || isreal(stanovenaKvalita)==0
        f = warndlg('Vektor "stanovenaKvalita" obsahuje nedovolené hodnoty (povolena pouze èísla 1, 2, 3).');
        return
    end
    
    podobnost=stanovenaKvalita;         
    podobnost(anotace~=stanovenaKvalita)=4;              % na pozice, kde se stanovená kvalita neshoduje s anotací se dá 4
    
    [jedna, dva, tri]=jdt(L,podobnost);
    TP=[jedna, dva, tri];            % poèet vzorkù v každé skupinì kvality
    clear pom jedna dva tri
    
    SE=TP./anotaceZastoupeni*100;% procentuální zastoupení správnì zaøazených vzorkù v každé skupinì
    
    [jedna, dva, tri]=jdt(L,stanovenaKvalita);
    kvalitaZastoupeni=[jedna, dva, tri];
    
    % výpoèet F1 skóre
    FP=kvalitaZastoupeni-TP;     % poèet falešnì pozitivních vzorkù v každé skupinì
    PP=TP./(TP+FP)*100; % pozitivní prediktivita
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
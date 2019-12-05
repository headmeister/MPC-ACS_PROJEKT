function [MAT]=principl_component(matice,p)
    % matice == vstupn� data (ve sloupc�ch)
    % p      == procentu�ln� zachov�n� informace
    % MAT    == v�stupn� matice
    
    
%% Standardizace sm�rodatou odchlkou
    % Z= [x - mean()]/std()
    % => nulov� pr�m�r a rozptyl roven 1
    matice=matice-mean(matice)./std(matice);
    
%% Normalizace
    % p�eps�n� do intervalu <-1;1>
    matice=matice./max(matice);

%% PCA
    C=cov(matice);
    [D,W]=eig(C);
    % D je matice vlastn�ch vektor� (v ��dc�ch)
    % W je matice na jej� uhlop���ce jsou vlastn� ��sla
    [rada(:,1),rada(:,2)]=sort(diag(W),'descend');
   %[hodnoty  ,jejich pozice]
    souc=sum(rada(:,1));
    rada_proc=rada(:,1)/souc;        % p�evedeno na procenta
    pod=1;
    while 1
        if sum(rada_proc(1:pod))>p   % kolik informace se p�enese do nov�ho syst�mu
            break
        else
            pod=pod+1;
        end
    end
    nov=rada(1:pod,:);
    
    MAT=matice*D(:,nov(:,2));         % p�eveden� do nov�ho sou�adnicov�ho syst�mu
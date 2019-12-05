function [MAT]=principl_component(matice,p)
    % matice == vstupní data (ve sloupcích)
    % p      == procentuální zachování informace
    % MAT    == výstupní matice
    
    
%% Standardizace smìrodatou odchlkou
    % Z= [x - mean()]/std()
    % => nulový prùmìr a rozptyl roven 1
    matice=matice-mean(matice)./std(matice);
    
%% Normalizace
    % pøepsání do intervalu <-1;1>
    matice=matice./max(matice);

%% PCA
    C=cov(matice);
    [D,W]=eig(C);
    % D je matice vlastních vektorù (v øádcích)
    % W je matice na jejíž uhlopøíèce jsou vlastní èísla
    [rada(:,1),rada(:,2)]=sort(diag(W),'descend');
   %[hodnoty  ,jejich pozice]
    souc=sum(rada(:,1));
    rada_proc=rada(:,1)/souc;        % pøevedeno na procenta
    pod=1;
    while 1
        if sum(rada_proc(1:pod))>p   % kolik informace se pøenese do nového systému
            break
        else
            pod=pod+1;
        end
    end
    nov=rada(1:pod,:);
    
    MAT=matice*D(:,nov(:,2));         % pøevedení do nového souøadnicového systému
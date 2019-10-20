function [vystup]=median_filtr(data,n,skok)
    % data   == vstupn� sign�l/matice, kter� projde filtrac� (ve sloupc�ch) 
    % n      == d�lka medi�nov�ho filtru
    % vystup == sign�l po filtraci
    vystup=medfilt1(data,n);
    vystup=sqrt(((vystup-data).^2)/length(data));
    vystup=resample(vystup,1,skok);

function [vystup]=median_filtr(data,n,skok)
    % data   == vstupní signál/matice, který projde filtrací (ve sloupcích) 
    % n      == délka mediánového filtru
    % vystup == signál po filtraci
    vystup=medfilt1(data,n);
    vystup=sqrt(((vystup-data).^2)/length(data));
    vystup=resample(vystup,1,skok);

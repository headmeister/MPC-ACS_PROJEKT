function[power]=spektrogram(data,skok) %vytvaøí výk. spektrogram bez pøekrývání z vektoru data po úsecích o délce skok

power=abs(spectrogram(data,skok,0,skok)).^2;

end
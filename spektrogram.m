function[power]=spektrogram(data,skok) %vytva�� v�k. spektrogram bez p�ekr�v�n� z vektoru data po �sec�ch o d�lce skok

power=abs(spectrogram(data,skok,0,skok)).^2;

end
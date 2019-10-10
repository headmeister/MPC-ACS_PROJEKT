function[rel_power]=drift_rel(power,skok,fvz,f_d,f_h) %vypo��t� relativn� mocnost frekven�n�ch slo�ek ze spektrogramu power,
%kde bloky maj� d�lku skok, vzorkovac� frekvence je fvz, f_d a f_h omezuj�
%rozsah jeho� relativn� mocnost chceme vypo��tat zad�v� se v Hz.
index=ceil(f_h/fvz*skok);
index2=ceil(f_d/fvz*skok);
n=1;
rel_power(1:length(power))=0;

while((n+skok)<length(power))
    rel_power(n:n+skok)=sum(power(n+index2:n+index))/sum(power(n:n+skok));
    n=n+skok;
end

index=ceil(f_h/fvz*(length(power)-n));
index2=ceil(f_d/fvz*(length(power)-n));

rel_power(n:end)=sum(power(n+index2:n+index))/sum(power(n:end));

end
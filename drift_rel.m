function[rel_power]=drift_rel(power,skok,fvz,f_d,f_h) %vypoèítá relativní mocnost frekvenèních složek ze spektrogramu power,
%kde bloky mají délku skok, vzorkovací frekvence je fvz, f_d a f_h omezují
%rozsah jehož relativní mocnost chceme vypoèítat zadává se v Hz.
index=ceil(f_h/fvz*skok);
index2=ceil(f_d/fvz*skok);

if(index2==0)
index2=1;
end
% n=1;
% rel_power(1:round(length(power)/skok))=0;
% ind=int32(1);
% while((n+skok)<length(power))
%     rel_power(ind)=sum(power(n+index2:n+index))/sum(power(n:n+skok-1));
%     n=n+skok;
%     ind=ind+1;
% end

rel_power=sum(power(index2:index,:))./sum(power);


% index=ceil(f_h/fvz*(length(power)-n));
% index2=ceil(f_d/fvz*(length(power)-n));
% 
% rel_power(ind)=sum(power(n+index2:n+index))/sum(power(n:end));

end
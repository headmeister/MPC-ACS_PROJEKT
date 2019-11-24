function[vyk]=accel(data,fvz_v,fvz_acc)



vyk=sqrt(data(1,:).^2+data(2,:).^2+data(3,:).^2);

data=data-mean(vyk);

vyk=sqrt(data(1,:).^2+data(2,:).^2+data(3,:).^2);
%   vyk=fft(vyk);
%   vyk(1:round(5/fvz_acc*length(data)))=0;
%   vyk(end+1-round(5/fvz_acc*length(data)):end)=0;
%   vyk=abs(ifft(vyk));
 %vyk=cumtrapz(vyk).^2;


 
vyk=diff(vyk);
filter=fir1(300,[2/(fvz_acc/5)],'low');

vyk=filtfilt(filter,1,vyk);
vyk=abs(hilbert(vyk));
 
%vyk=vyk.^(1/4);
vyk=resample(vyk,fvz_v,fvz_acc);



end
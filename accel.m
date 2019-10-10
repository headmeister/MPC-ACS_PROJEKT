function[vyk]=accel(data,fvz_v,fvz_acc)



vyk=sqrt(data(1,:).^2+data(2,:).^2+data(3,:).^2);
vyk=resample(vyk,fvz_v,fvz_acc);
vyk=fft(vyk);
vyk(1)=0;
vyk=ifft(vyk);


end
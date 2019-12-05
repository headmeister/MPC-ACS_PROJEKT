clear all
P1=load('P1.mat');
P2=load('P2.mat');
P3=load('P3.mat');
P12=load('P12.mat');
P23=load('P23.mat');
P123=load('P123.mat');
ch(1:9)=0;
progress=0;
fis=readfis('FIS_std_fce.fis');

an3=load('111001_ann_quality_train.mat');
an2=load('105001_ann_quality_train.mat');
an1=load('100001_ann_quality_train.mat');

anotace=[an1.anotace an2.anotace an3.anotace];
anotace_pom=round(resample(anotace,1,1000));
for i=1:9
 ch(i)=mse(P3.P(i,:)-P123.P(i,end-length(P3)+1:end));   



end


AN1=round(resample(an1.anotace,1,1000));
AN1=AN1(1:length(P1.P));
AN2=round(resample(an2.anotace,1,1000));
AN2=AN2(1:length(P2.P));
AN3=round(resample(an3.anotace,1,1000));
AN3=AN3(1:length(P3.P));


P=[P1.P P2.P P3.P];
AN=[AN1 AN2 AN3];
AN(72000:95000)=[];
AN=ind2vec(AN);
PP=P;
PP(:,72000:95000)=[];

%% trenovani
sit=patternnet([12 12]); % vytvoøení sítì s n skrytými neurony
sit.trainParam.max_fail=500;
sit.trainParam.mu=0.0001;
sit.trainParam.mu_dec=0.0001;
sit.trainParam.mu_inc=0.0002;
sit.trainParam.epochs=10000;
sit=train(sit,PP,AN,'useParallel','yes'); % uèení sítì na trénovací množinì
FF=sit(P,'useParallel','yes'); % vybavení sítì ze vstupu pøíznakù
%% pokracovani
    kval=FF(1,:)+FF(2,:)*2+FF(3,:)*3;
    %kval=[1 kval 1];
disp('vybaveno'); 
    
    
  
   mat(1:length(kval)-2,1:3)=0;
  for(i=2:length(kval)-1)
      mat(i-1,:)=kval(i-1:i+1);
      
      
  end  
  
kval(2:end-1)=evalfis(mat,fis);
  
  
%kval=[pom(1)', kval', pom(end)'];
    
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
disp('FIS evaluace kompletní');    
 %kval= vec2ind(round(FF)); % návrat ke kvalitì
 kval_v=(resample(kval,round(length(anotace)/length(FF(1,:))),1)); % pøevzorkování zpìt
% 
%kval=medfilt1(round(kval),3*skok); % spojení do delších úsekù
kval_v(kval_v>3)=3;
kval_v(kval_v<1)=1;


kval_v(kval_v==0)=1;
kval_v=round(kval_v);

if(length(kval_v)>length(anotace))
    ppp=length(kval_v)-length(anotace);
    kval_v=kval_v(1:end-ppp);
elseif (length(kval_v)<length(anotace))
    ppp=length(anotace)-length(kval_v);
    ppp=kval_v(end)*ones(1,ppp);
    kval_v=[kval_v ppp];
end
progress=progress+10;
disp(['Progress: ' num2str(progress) '%']);
disp('Hotovo');   
%  [uspesnost, uspescnost_cn]=hodnoceni(kval_v,anotace);

[SE, PP, F1, cov]=F1score(anotace,kval_v);
disp(['SE: ' num2str(SE,3) ' %' ]);
disp(['PP: ' num2str(PP,3) ' %']);
disp(['F1: ' num2str(F1,3) ', coverage: ' num2str(cov,3) ' %']);



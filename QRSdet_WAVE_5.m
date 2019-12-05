function QRS = QRSdet_WAVE_5(fvz,vstup_sig)
%------ NAPOVEDA ----------------------------------------------------------

% vstupni parametry
% fvz : pouzita vzorkovaci frekvence v Hz
% vstup_sig : vstupni signal ve formatu *.mat, radkovy vektor

% vystupni parametry
% QRS : pozice QRS komplexu v zadanem svodu

%------ NASTAVENI PARAMETRU -----------------------------------------------

% parametry programu
typ=2;

%         k(1)  % pásmo
%         k(2)  % násobek prahu	
%         k(3)  % maximalni povolena vzdalenost mezi dvojici extremu v sekundach
%         k(4)  % minimalni vzdalenost mezi dvema QRS komplexy v sekundach
%         k(5)  % délka okna pro výpoèet prahu
%         k(6)  % násobek druhého prahu
%         k(7)  % velikost okraje
%         k(8)  % maximální vzdálenost dvou QRS

if typ == 1
   k = [28 1.803 0.05 0.4 1.91*fvz 1.19 0.45 1.2*fvz];  % 98.9550%
elseif typ == 2
   k = [15.1 1.73 0.112 0.29 2.05*fvz 1.2 0 1.12*fvz];  % 99.6038%
else
   k = [22 1.85 0.112 0.29 2.05*fvz 1.2 0 1.12*fvz];    % 99.0159%
end

pasmo = k(1);
vlnka = 'bior1.5';      % typ vlnky

%------ CWT ---------------------------------------------------------------

band = transformace(vstup_sig,pasmo,vlnka);

%------ VYPOCET PRAHU -----------------------------------------------------

prah = movingstd(band,round(k(5)/2),'central');

%------ DETEKCE QRS -------------------------------------------------------

QRS = polohy_vln(band,k(2)*prah,fvz*k(3)); % funkce nalezne polohy QRS v pasmu

%------ ELIMINACE MULTIVLN ------------------------------------------------

elim = diff(QRS);                   
kde = find(elim < fvz*k(4));
QRS(kde+1) = [];

%------ KOREKCNI PRAVIDLA -------------------------------------------------

again = find(diff(QRS)>k(8));
kraj = round(k(7)*fvz);

new = [];
for i = 1:length(again)
    
    od = QRS(again(i))+kraj;
    do = QRS(again(i)+1)-kraj;
    
    pom = polohy_vln(band(od:do),k(6)*prah(od:do),fvz*k(3));
    new = [new pom+od-1];
    
end

QRS = [QRS new];
QRS = sort(QRS);

elim = diff(QRS);                   
kde = find(elim < fvz*k(4));
QRS(kde+1) = [];

end
function [bod] = polohy_vln(x,prah,mez)
%----- INICIALIZACE -------------------------------------------------------
  
bod = [];
extrem_m = [];
extrem_p = [];
   
%------ NALEZENÍ PRÙCHODÙ NULOU -------------------------------------------

u = find(diff(sign(x)));
u = [1 u length(x)];

%------ NALEZENÍ EXTREMU --------------------------------------------------

x1 = find(x>=prah);
poz1 =  find(diff(x1)-1);
poz1 = [0 poz1 length(x1)];
if length(poz1) > 2
    for i = 1:length(poz1)-1
        k = x( x1(poz1(i)+1:poz1(i+1)) );
        extrem_p(i) = find(k == max(k)) + x1(poz1(i)+1)-1;
    end
end

x2 = find(x<=-prah);
poz2 =  find(diff(x2)-1);
poz2 = [0 poz2 length(x2)];
if length(poz2) > 2
    for i = 1:length(poz2)-1
        k = x( x2(poz2(i)+1:poz2(i+1)) );
        extrem_m(i) = find(k == min(k)) + x2(poz2(i)+1)-1;
    end
end

%------ NALEZENÍ PRUCHODU NULOU MEZI EXTREMY ------------------------------

if ~isempty(extrem_p) && ~isempty(extrem_m)
    
for i = 1:length(extrem_p)
    k = find(extrem_m>extrem_p(i),1,'first');
    if ~isempty(k)
        kk = find(extrem_p<extrem_m(k),1,'last');
        if  extrem_p(i)==extrem_p(kk) && extrem_m(k)-extrem_p(i)< mez
            poloha = find(u>=extrem_p(i) & u<=extrem_m(k));
            poloha = u(poloha(ceil(length(poloha)/2)));
        bod = [bod poloha];
        end
    end
end

for i = 1:length(extrem_m)
    k = find(extrem_p>extrem_m(i),1,'first');
    if ~isempty(k)
        kk = find(extrem_m<extrem_p(k),1,'last');
        if  extrem_m(i)==extrem_m(kk) && extrem_p(k)-extrem_m(i)< mez
            poloha = find(u>=extrem_m(i) & u<=extrem_p(k));
            poloha = u(poloha(ceil(length(poloha)/2)));
        bod = [bod poloha];
        end
    end
end

end
%------ SERAZENI POZIC ----------------------------------------------------

bod = sort(bod);

end
function [band] = transformace(vstup_sig,pasmo,vlnka)
%------ NAPOVEDA ----------------------------------------------------------

% funkce provadi rychlou CWT s vyuzitim konvoluce

%------ VYTVORENI FILTRU --------------------------------------------------

if  pasmo == 15 && strcmp(vlnka,'bior1.5')
    imp_char = [-7.60565599708887e-16,-1.85659367047585e-11,-4.45140058022962e-09,8.97251273307321e-08,-6.53831413699163e-07,-3.58815479937255e-06,1.33324726230374e-05,1.72677727375859e-05,-8.13034241591086e-05,-0.000146486039457173,-0.000158209850004651,0.000145613639643926,0.00107927082419663,0.00121693755225996,0.000132862465243511,-0.00200069372513942,-0.00517203607827470,-0.00658411822754586,-0.00590073070315290,-0.00683068849394393,-0.00668885498518011,0.000415044740183573,0.0138419336827928,0.0364237444725993,0.0550911264396067,0.0562097611661096,0.0464790430912854,0.0309942198664056,-0.00410399118501448,-0.0576403173228762,-0.135672219588872,-0.238644260830236,-0.303562256437038,-0.331380919672457,-0.340130389822616,-0.308264798217232,-0.247467189706591,-0.152279038815097,-0.00840605215956051,0.139800722880035,0.239634586842784,0.303380284390641,0.338943983730319,0.332878660953073,0.307990683019033,0.244958260784066,0.148316651348171,0.0643972223407127,0.00884067750687203,-0.0285538987166462,-0.0452044560796249,-0.0554507879029005,-0.0560798652309547,-0.0389972661784107,-0.0156749561683949,-0.00145430879176865,0.00624697761669970,0.00702698923766024,0.00588666001736941,0.00655645951063974,0.00538761083014166,0.00235368499262149,3.87916588733341e-05,-0.00115460182689205,-0.00116549674890849,-0.000213031569019047,0.000152489701845874,0.000145536489723662,9.46800031922444e-05,-1.18267930921447e-05,-1.69189414038318e-05,3.42459175725230e-06,1.09439007460723e-06,-1.76864434296067e-07,1.25969867651953e-08,2.08973938848534e-10];
elseif pasmo == 41 && strcmp(vlnka,'bior1.5')    
    imp_char = [-1.25742756689685e-15,-1.26881392334892e-13,-8.73922309057134e-12,-2.88645843074840e-10,2.97511947557746e-10,5.26727722869422e-09,-3.85925438746758e-08,-4.98158678247319e-08,1.38768100942093e-07,4.17137830605845e-07,-1.05643472903079e-07,-1.38865815307360e-06,-2.13690043263204e-06,-2.09377032496627e-06,-1.97023584679570e-06,2.39748851998860e-06,1.23486742999098e-05,1.94488878778454e-05,1.51784277959329e-05,2.52125057878260e-06,-1.79340650284360e-05,-4.84485880654922e-05,-7.97602833560836e-05,-9.47492497252575e-05,-8.66161406134489e-05,-8.45116725592607e-05,-9.83820978497574e-05,-0.000105001790035676,-4.77225496081356e-05,4.82895244543950e-05,0.000204124685821793,0.000433427530170947,0.000699583234529610,0.000834758415229956,0.000826072378541680,0.000712739837390467,0.000523355623062375,0.000238704928453196,-0.000142000883153561,-0.000582574208477941,-0.00111678631390131,-0.00170671966780629,-0.00255294069225966,-0.00321453155356668,-0.00371366380068779,-0.00402910409422038,-0.00406701277460427,-0.00377511700943481,-0.00356401342128039,-0.00355623518565073,-0.00367967556443057,-0.00398034003646380,-0.00454347381780396,-0.00474231981668435,-0.00411038585644414,-0.00270336142960599,-0.000953120375467694,0.00115742563169297,0.00386087844623606,0.00730828699222139,0.0112316524688249,0.0160393029098479,0.0218181037465998,0.0275732902310478,0.0317314832059487,0.0344470636794224,0.0356636029416735,0.0351534851762619,0.0332361373995689,0.0312156560190402,0.0288557167904198,0.0260202584551251,0.0227187280935774,0.0187504961487219,0.0126907456233335,0.00427873444741241,-0.00573890534619091,-0.0167496344934771,-0.0290146064179105,-0.0438370354098923,-0.0606339181897053,-0.0796755004643571,-0.0969583674084762,-0.124726758030298,-0.146417967009399,-0.164232250289230,-0.178378314311040,-0.189240541420099,-0.195938838583113,-0.200515494199318,-0.203806466540297,-0.205403214959108,-0.204504484593999,-0.201416165839947,-0.195558086652470,-0.186386706649258,-0.173906028449381,-0.159464825486535,-0.143307411178697,-0.123984058705305,-0.100204426861308,-0.0730368340364784,-0.0415330650695232,-0.00537037651685421,0.0319251040679111,0.0646461208167127,0.0929789916496002,0.117701454228338,0.138363090846845,0.155101757772216,0.170043987858716,0.183166383615305,0.193342817648000,0.200095236739690,0.203853723034840,0.205403880218856,0.204460255706795,0.201587133006526,0.197317641128100,0.191524183895262,0.181784529058429,0.168540447735490,0.151769347058811,0.126448539430783,0.108151635804853,0.0854955007483108,0.0656785442975078,0.0483923171104386,0.0328834482711410,0.0200122447567570,0.00871412569850927,-0.00156949938283585,-0.0105801216638661,-0.0172817797255345,-0.0217481857604045,-0.0251200124310823,-0.0281250518503231,-0.0305923260386450,-0.0326849045749701,-0.0346764028796881,-0.0357244909098922,-0.0349314750193850,-0.0326458100138954,-0.0288894088423525,-0.0235464705580775,-0.0175487556090831,-0.0124859514843156,-0.00834430845645248,-0.00476978783852583,-0.00183243434276977,0.000400430378799387,0.00224235561300576,0.00377814994283247,0.00464850489835909,0.00468088476109749,0.00411175320161333,0.00374607536966372,0.00357149738534964,0.00355166657327058,0.00368394086355591,0.00401466535109187,0.00406948640488193,0.00382800248803475,0.00336561452306181,0.00266532787104792,0.00200986318495701,0.00129236369301020,0.000717493087728296,0.000261858890739376,-0.000141753035798926,-0.000452970416298559,-0.000669020848803373,-0.000801000985317976,-0.000847059653789262,-0.000752557340832607,-0.000512580097557159,-0.000259523835261330,-8.50711939328661e-05,2.44545414162085e-05,9.48245104263474e-05,0.000104460291201567,8.64898764976676e-05,8.42755622739342e-05,9.43739310153585e-05,8.59290885945942e-05,5.83333604731448e-05,2.51867213902288e-05,2.53599482507695e-06,-1.25248535268426e-05,-1.91860853303255e-05,-1.54475077203722e-05,-4.54267790716980e-06,1.20319401187173e-06,2.27685253125787e-06,2.07732346221462e-06,1.76458893852726e-06,3.92036692800090e-07,-3.51672021222431e-07,-2.60408465211960e-07,3.64562732130277e-08,4.61779802462585e-08,2.74283112596608e-09,-3.35123314764449e-09,6.55357819646668e-10,-4.61460502880443e-11,-8.09747680208346e-13];
else
    impulz = [zeros(1,round(pasmo*5)) 1 zeros(1,round(pasmo*5))];
    odezva = cwt(impulz,pasmo,vlnka);
    imp_char = odezva(find(odezva));
end
%------ KONVOLUCE ---------------------------------------------------------

zp = (length(imp_char)-1)/2;
zac(1:ceil(zp)) = vstup_sig(1);
kon(1:ceil(zp)) = vstup_sig(end);

vstup_sig = [zac vstup_sig kon];

band = conv(vstup_sig,imp_char);

%------ ELIMINACE ZPOZDENI ------------------------------------------------

band = band(ceil(zp)+1:end-fix(zp));
band = band(ceil(zp)+1:end-ceil(zp));

end
function s = movingstd(x,k,windowmode)
% movingstd: efficient windowed standard deviation of a time series
% usage: s = movingstd(x,k,windowmode)
%
% Movingstd uses filter to compute the standard deviation, using
% the trick of std = sqrt((sum(x.^2) - n*xbar.^2)/(n-1)).
% Beware that this formula can suffer from numerical problems for
% data which is large in magnitude. Your data is automatically
% centered and scaled to alleviate these problems.
%
% At the ends of the series, when filter would generate spurious
% results otherwise, the standard deviations are corrected by
% the use of shorter window lengths.
%
% arguments: (input)
%  x   - vector containing time series data
%
%  k   - size of the moving window to use (see windowmode)
%        All windowmodes adjust the window width near the ends of
%        the series as necessary.
%
%        k must be an integer, at least 1 for a 'central' window,
%        and at least 2 for 'forward' or 'backward'
%
%  windowmode - (OPTIONAL) flag, denotes the type of moving window used
%        DEFAULT: 'central'
%
%        windowmode = 'central' --> use a sliding window centered on
%            each point in the series. The window will have total width
%            of 2*k+1 points, thus k points on each side.
%        
%        windowmode = 'backward' --> use a sliding window that uses the
%            current point and looks back over a total of k points.
%        
%        windowmode = 'forward' --> use a sliding window that uses the
%            current point and looks forward over a total of k points.
%
%        Any simple contraction of the above options is valid, even
%        as short as a single character 'c', 'b', or 'f'. Case is
%        ignored.
%
% arguments: (output)
%  s   - vector containing the windowed standard deviation.
%        length(s) == length(x)

% check for a windowmode
if (nargin<3) || isempty(windowmode)
  % supply the default:
  windowmode = 'central';
elseif ~ischar(windowmode)
  error 'If supplied, windowmode must be a character flag.'
end
% check for a valid shortening.
valid = {'central' 'forward' 'backward'};
ind = find(strncmpi(windowmode,valid,length(windowmode)));
if isempty(ind)
  error 'Windowmode must be a character flag, matching the allowed modes: ''c'', ''b'', or ''f''.'
else
  windowmode = valid{ind};
end

% length of the time series
n = length(x);

% check for valid k
if (nargin<2) || isempty(k) || (rem(k,1)~=0)
  error 'k was not provided or not an integer.'
end
switch windowmode
  case 'central'
    if k<1
      error 'k must be at least 1 for windowmode = ''central''.'
    end
    if n<(2*k+1)
      error 'k is too large for this short of a series and this windowmode.'
    end
  otherwise
    if k<2
      error 'k must be at least 2 for windowmode = ''forward'' or ''backward''.'
    end
    if (n<k)
      error 'k is too large for this short of a series.'
    end
end

% Improve the numerical analysis by subtracting off the series mean
% this has no effect on the standard deviation.
x = x - mean(x);
% scale the data to have unit variance too. will put that
% scale factor back into the result at the end
xstd = std(x);
x = x./xstd;

% we will need the squared elements 
x2 = x.^2;

% split into the three windowmode cases for simplicity
A = 1;
switch windowmode
  case 'central'
    B = ones(1,2*k+1);
    s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/(2*k+1)))/(2*k));
    s(k:(n-k)) = s((2*k):end);
  case 'forward'
    B = ones(1,k);
    s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/k))/(k-1));
    s(1:(n-k+1)) = s(k:end);
  case 'backward'
    B = ones(1,k);
    s = sqrt((filter(B,A,x2) - (filter(B,A,x).^2)*(1/k))/(k-1));
end

% special case the ends as appropriate
switch windowmode
  case 'central'
    % repairs are needed at both ends
    for i = 1:k
      s(i) = std(x(1:(k+i)));
      s(n-k+i) = std(x((n-2*k+i):n));
    end
  case 'forward'
    % the last k elements must be repaired
    for i = (k-1):-1:1
      s(n-i+1) = std(x((n-i+1):n));
    end
  case 'backward'
    % the first k elements must be repaired
    for i = 1:(k-1)
      s(i) = std(x(1:i));
    end
end

% catch any complex std elements due to numerical precision issues.
% anything that came out with a non-zero imaginary part is
% indistinguishable from zero, so make it so.
s(imag(s) ~= 0) = 0;

% restore the scale factor that was used before to normalize the data
s = s.*xstd;
end
t=Demeter(:,1); 
u=Demeter(:,2);
y=Demeter(:,3);
x=Demeter(:,4);

figure
plot(t,u,t,y)
legend('u','y')

%%
i1=379; %intrare(u)
i2=387;
i3=381; %iesire(y)
i4=390;
K = mean(y)/mean(u) %raportul dintre media ieșirii și media intrării
%i5=59;
%i6=91;
%i7=61;
%i8=94;
Mr = (y(i4)-y(i3))/(u(i2)-u(i1)/K) %Amplificare maxima
%tita = (sqrt(1-(sqrt(Mr^2*(Mr^2-1))/Mr^2)))/sqrt(2); 
tita=sqrt((Mr-sqrt(Mr^2-1))/(2*Mr)) %factorul de amortizare
Tr=2*(t(i2)-t(i1)) %Perioada de rezonanta 
wr=2*pi/Tr  %Pulsatia de rezonanta
wn=wr/sqrt(1-2*tita^2) %Pulsatia naturala
 
H= tf(K*wn^2,[1 2*tita*wn wn^2]) %functia de transfer

%validarea datelor folosind spatiul starilor
A=[0,1;-wn^2,-2*tita*wn];
B=[0;K*wn^2];
C=[1,0];
D=0;
 
sys=ss(A,B,C,D);
ysim=lsim(sys,u,t,[y(1),(y(2)-y(1))/(t(2)-t(1))]);


J=norm(y-ysim)/sqrt(length(y)) %eroarea medie patratica
Empn=norm(y-ysim)/norm(y-mean(y))*100 %eroarea medie patratica normalizata

figure
plot(t,ysim,t,y)
legend('Ieșirea simulată','y')
hold on

%% identificare armax
%Te=t(2)-t(1);
%data_id = iddata(y(i5:i6),u(i5:i6),Te);
m_armax = armax(data_id,[2 2 2 1]);

%gradul de suprapunere
figure;compare(data_id,m_armax);shg
title("Metoda armax compare")
% validarea statistica
figure;resid(data_id,m_armax)
title("Metoda armax resid")
%functia de transfer
Ha=tf(m_armax.B,m_armax.A,Te,'variable','z^-1');
zpk(Ha)
%%
%perioada de achizitie/esantionare
%data_id=iddata(y(i5:i6),u(i5:i6),Tr);
%data_vd=iddata(y(i7:i8),u(i7:i8),Tr);
%% OE

nF=2;%un pol
nB=2;%nu avem zerouri
nD=1;%un tact de intarziere
m_oe=oe(data_id,[nB,nF,nD])

%validare statistica
figure; 
resid(m_oe,data_vd)

%gradul de suprapunere
figure; compare(m_oe,data_vd)
dt = t(2)-t(1)
H4 = tf(K*wn^2,[1 +2*tita*wn wn^2])
c2d(H4,dt,'zoh')


%% identicare cu iv
%perioada de achizitie/esantionare

data_id=iddata(y(i5:i6),u(i5:i6),Tr);
data_vd=iddata(y(i7:i8),u(i7:i8),Tr);

%metoda iv
nA=2;
nB=1;
nd=1;
m_iv=iv4(data_id,[nA,nB,nd]);

%val statisctica
figure
resid(m_iv,data_vd);

%gradul de suprapunere
figure
compare(m_iv,data_vd);




close all;clear all;clc;
%stale 
c = 2.99792458*1e8; mu0 = 4 * pi * 1e-7; eps0 = 1/(mu0*c*c);

% parametry symulacj
a=0.02;
d=0.4;
M=10;
N=200;
x=linspace(0,a,M)';
dx=x(2)-x(1);
y=linspace(0,d,N)';
dy=y(2)-y(1);
I=M*N;

% wpolrzedne wezlow siatki
[X,Y]=meshgrid(x,y); X=X'; Y=Y';

% warunki poczatkowe
e0=zeros(I,1);
hx0=zeros(I,1);
hy0=zeros(I,1);

% warunki brzegowe
temp=zeros(M,N);
temp(1,:)=1;
temp(M,:)=1;
temp(:,1)=1;
temp(:,N)=1;
bnde=find(temp==1);
%bnde=[bnde ;803];

temp=zeros(M,N);
temp(1,:)=1;
temp(M,:)=1;
temp(:,N)=1;
bndhx=find(temp==1);

temp=zeros(M,N);
temp(M,:)=1;
temp(:,1)=1;
temp(:,N)=1;
bndhy=find(temp==1);

% okreslenie czasu analizy
dt=dx*dy/c/(dx+dy)/2; % krok czasowy
K=10000;  % zatem calkowity czas analizy = dt*K

% pobudzenie 
f=8e9; 
w=2*pi*f;
t=dt*(1:K);
J_0=1;
tc=dt*K/25;
ts=dt*K/50;
J_exc=J_0*sin(t*w).*exp(-((t-tc)/ts).^2); % impulsowe 
J_exc(t>tc)=J_0*sin(t(t>tc)*w); % przechodzi w ciagle po odkomentowaniu
m_e=3; % lokalizacja pobudzenia
n_e=10;
i_e=n_e*M+m_e;
 
% warunek radiadcyjny
beta=pi/a;
v=w/sqrt((w/c)^2-beta^2);
s=(v*dt-dy)/(v*dt+dy);
 
% obliczenia - iteracje
hnde=[]; figure(1); hold on;
for k=1:K

hx = hx0 - dt / mu0 / dy * ( [e0(M+1:I);zeros(M,1)] - e0 );
hy = hy0 + dt / mu0 / dx * ( [e0(2:I);0] - e0 );
hx(bndhx)=0;
hy(bndhy)=0;

e = e0 + dt/eps0 * ( ( hy - [0;hy(1:I-1)] )/dx - (hx - [zeros(M,1);hx(1:I-M)])/dy);
e(bnde)=0; % PEC wszedzie 
e(I-M+1:I)= e0(I-2*M+1:I-M)+s*( e(I-2*M+1:I-M)-e0(I-M+1:I) ); % otwarty koniec

e(i_e)=e(i_e)-dt/eps*J_exc(k);%
delete(hnde); hnde=mesh(X,Y,reshape(e,M,N)); xlim([0 d]); ylim([0 d]); zlim([-1 1]*2e4); drawnow;
e0=e;
hx0=hx;
hy0=hy;
end



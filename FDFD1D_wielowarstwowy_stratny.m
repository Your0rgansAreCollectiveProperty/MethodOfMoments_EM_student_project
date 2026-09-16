clear; close all; clc;
%stale
c = 2.99792458*1e8;
mu0 = 4 * pi * 1e-7;
eps0 = 1/(mu0*c*c);


% czestotliwosc
f=(183+0i)*1e6;%c/L/2; 
w=2*pi*f;

% parametry symulacji
L=1;
M=1000;
z=linspace(0,L,M);
dz=z(2)-z(1); % krok przestrzenny

eps=ones(M,1);
mu=ones(M,1);
 
eps((z>0.7)&(z<0.9))=4;
mu((z>0.7)&(z<0.9))=4;

sigma=1;
eps((z>0.9))=1+sigma/(1j*w*eps0);

eps=eps0*diag(eps);
mu=mu0*diag(mu);





% konstrukcja operatorow macierzowych
RE=(diag(-ones(M,1),0)+diag(ones(M-1,1),1))/dz;
RH=RE';
A=RH*mu^(-1)*RE-w^2*eps*eye(M,M);

% pobudzenie
J_exc=1; % ciagle 
dist=0.1;% lokalizacja pobudzenia
exc_in=round(dist/dz); % wezel w ktorym zlokalizowane jest pobudzenie
B=zeros(M,1);
B(exc_in,1)=-1j*w*J_exc;

% warunki brzegowe
A([1,M],:)=[];
A(:,[1,M])=[];
B([1,M],:)=[];

% rozwiazanie zagadnienia
e=A\B;
e=[0;e;0];


% wizualizacja
Emax=max(abs(e))
T=1/f;
for t=0:T/100:10*T
    plot(z,real(e*exp(1j*w*t))); 
    ylim([-1 1]*Emax);
    drawnow;    
end


%============================================
% Curso de Comunicaciones Opticas Coherente
% Dr. Mario R. Hueda
% Contenido: Medicion del Vector de Stokes
% Version: 24_10_17
%============================================
close all
clear

OS=200; 
wt=[-1000*OS:1:1000*OS]/OS*pi;

% Angulos del Vector de Jones ("deterministico")
phi=pi/4*1.0;
chi=pi/3*1.0;

% Angulos del Vector de Jones ("aleatorio")
%phi=2*rand*pi;
%chi=rand*pi;

phi=mod(phi,2*pi)-pi;
chi=mod(chi,2*pi)-pi;

%--------------------------------------------------------------------
% Vector de Stokes "teorico"
%--------------------------------------------------------------------
S_teo = [cos(2*chi); sin(2*chi)*cos(phi); sin(2*chi)*sin(phi)];

%--------------------------------------------------------------------
% Componentes temporales del campo electrico a partir del vector de Jones
%--------------------------------------------------------------------
a=exp(j*pi*.0); % cambia el vector de Stokes al varia "a"?
E_x=a*cos(chi)*cos(wt);
E_y=a*sin(chi)*cos(wt+phi);

%--------------------------------------------------------------------
% Medicion de intensidades en "x" e "y"
%--------------------------------------------------------------------
Ix=mean(E_x.*conj(E_x));
Iy=mean(E_y.*conj(E_y));

%--------------------------------------------------------------------
% Medicion de intensidades a "+45" y "-45"
%--------------------------------------------------------------------
tita=1.0*pi/4;
e1=E_x.*cos(tita)^2+E_y.*sin(tita)*cos(tita);
e2=E_y.*sin(tita)^2+E_x.*sin(tita)*cos(tita);
I_45=real(mean(conj(E_x).*e1+conj(E_y).*e2));

tita=-1.0*pi/4;
e1=E_x.*cos(tita)^2+E_y.*sin(tita)*cos(tita);
e2=E_y.*sin(tita)^2+E_x.*sin(tita)*cos(tita);
I_m45=real(mean(conj(E_x).*e1+conj(E_y).*e2));

%--------------------------------------------------------------------
% Medicion de intensidades en waveplate de lambda/4 rotado a 45
%--------------------------------------------------------------------

% a) Rotacion de ejes (alineacion con ejes ppales del "waveplate")
alfa=1.0*pi/4;
E_xr=E_x.*cos(alfa)+E_y.*sin(alfa);
E_yr=-E_x.*sin(alfa)+E_y.*cos(alfa);

% b) Desfasaje en pi/2 para una componente (linea de retardo)
E_xr=[zeros(1,OS/2) E_xr(1:end-OS/2)];

% c) "Des-rotacion" de ejes (alineacion con ejes de coordenadas)
alfa=-1.0*pi/4;
e1=E_xr.*cos(alfa)+E_yr.*sin(alfa);
e2=-E_xr.*sin(alfa)+E_yr.*cos(alfa);
E_xr=e1(OS:end);% elimina primeras OS muestras p/mejorar estimacion
E_yr=e2(OS:end);% elimina primeras OS muestras p/mejorar estimacion

% Medicion de intensidades en "x" e "y"
I_R=real(mean(conj(E_xr).*E_xr));
I_L=real(mean(conj(E_yr).*E_yr));

%--------------------------------------------------------------------
% "Medicion" del Vector de Stokes
%--------------------------------------------------------------------
s0=Ix+Iy;
s1=Ix-Iy;
s2=I_45-I_m45;
s3=I_R-I_L;
S_med=[s1 s2 s3]'/s0;

%--------------------------------------------------------------------
% "Medicion" del Vector de Jones
%--------------------------------------------------------------------
s_dot_sigma=[S_med(1), S_med(2)-j*S_med(3);S_med(2)+j*S_med(3), -S_med(1)];
[V v]=eig(s_dot_sigma);
m=find(abs(diag(v)-1)<1e-5); % Busqueda del autovector corresp. al autovalor=1
J=V(:,m);
J=exp(-j*angle(J(1)))*J; % Vector Jones

% Angulos del vector de Jones
phi_med=angle(J(2));
chi_med=atan(abs(J(2))/real(J(1)));

%--------------------------------------------------------------------
% Vector de Stokes a partir de Jones con los angulos medidos
%--------------------------------------------------------------------
S_med2 = [cos(2*chi_med); sin(2*chi_med)*cos(phi_med); sin(2*chi_med)*sin(phi_med)];

%--------------------------------------------------------------------
% Presentacion de resultados
%--------------------------------------------------------------------
fprintf('\nVector de Stokes (teoria)=[%1.2f, %1.2f, %1.2f]',S_teo(1),S_teo(2),S_teo(3));
fprintf('\nVector de Stokes (medido)=[%1.2f, %1.2f, %1.2f]',S_med(1),S_med(2),S_med(3));
fprintf('\nVector de Stokes (medido)=[%1.2f, %1.2f, %1.2f]\n',S_med2(1),S_med2(2),S_med2(3));


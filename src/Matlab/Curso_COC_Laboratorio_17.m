%=======================================================================
% Comunicaciones Opticas Coherentes
% Prof.: Dr. Mario Hueda (mario.hueda@unc.edu.ar)
% Practico Laboratorio
%=======================================================================
clear all;
close all

%============================================
% Ejemplo 1: para PLL en baud-rate TR (e.g., M&M): BR=42 - D=17 - Bl=1.25e6 - P=160
% Ejemplo 2: para Low Latency PLL:  BR=42 - D=1 - Bl=50e6 - P=40 - Ki=0
%============================================

% Frecuencia de simbolos [Hz]
BR=42.0e9;

% Frecuencia real de trabajo [Hz]
P=40; % Factor de paralelismo
fs=BR/P; 

% Ancho de banda del lazo [Hz]
Bl=50.e6;

% Ganancias del lazo
Kp=2*pi*Bl/fs;   % Proporcional
Ki=Kp/1000*.0; % Integral

% Latencia del lazo en ciclos de reloj
D=1;  


%============================================
% Evaluacion de la Funcion de Transferencia
%============================================
N=10000;                 % Numero de puntos de la respuesta en frecuencia entre 0 y fs/2
w=logspace(-10,pi,N+1);
f=(0.5*fs/pi)*w; % Nota: w=2*pi*f*fs
z=exp(j*w);

% Respuestas de lazo abierto y cerrado
Lz=z.^(-D).*(Kp+Ki./(1-z.^(-1)));
Gz=Lz.*z.^(-1)./(1-z.^(-1));
Hz=Gz./(1+Gz);
Gdb=20.0*log10(abs(Hz));
phase=unwrap(angle(Hz));

figure(1)
subplot 211
h=semilogx(f/1e6,Gdb,'r',f/1e6,-3*ones(1,length(f)),'-b',f/1e6,1*ones(1,length(f)),'-b');
set(h,'Linewidth',3);
set(gca,'XScale','log','YScale','lin','FontWeight','bold','FontSize',18);
set(gca,'Linewidth',2);
axis([1.0e-2,200,-20,4]);
xlabel('Frequency [MHz]');
ylabel('Amplitude [dB]');
title('Jitter Transfer Function');
grid on;

subplot 212
h=semilogx(f/1e6,phase,'r');
set(h,'Linewidth',3);
set(gca,'XScale','log','YScale','lin','FontWeight','bold','FontSize',18);
set(gca,'Linewidth',2);
axis([1.0e-2,200,-4,1]);
xlabel('Frequency [MHz]');
ylabel('Phase [rad]');
grid on;

%=======================================================================
% Comunicaciones Opticas Coherentes
% Prof.: Dr. Mario Hueda (mario.hueda@unc.edu.ar)
% Practicos Lab. 2-5
%=======================================================================
clear;
close all;
ON=1;
OFF=0;

%============================================
% Generacion de la Respuesta al Impulso
%============================================
fB = 32e9;	% Velocidad de simbolos (baud rate)
T = 1/fB; % Tiempo entre simbolos
M = 4;  %Factor de sobremuestreo
fs = fB*M;	% Sample rate
ModOrd=64; %"4"=QPSK / "16"=QAM16

SNRdB=25;
SNR=10^(SNRdB/10);
n_symbols = 500000;
flag_mod_sim_2=ON; %"Off" Modelo 1 / "ON" Modelo 2

%============================================
% Respuesta del Filtro Transmisor
%============================================
beta = 0.1; %Factor de roll-off
gn=rcosine(fB,fs,'sqrt',beta,24);  %  Generacion usando funcion de matlab (raiz cuadrada)
gn2=rcosine(fB,fs,'normal',beta,24);  %  Generacion usando funcion de matlab (coseno realzado)

%============================================
% Respuesta del Filtro Transmisor+Canal (h[n]=b[n]*g[n])
%============================================
bn=1;
hn=conv(gn,bn);
hn=hn/sum(hn)*M;
if flag_mod_sim_2==1
    hn=hn/sqrt(sum(hn.^2));
end
Eh=sum(hn.^2);

%============================================
% Respuesta del Filtro Receptor (f[n])
%============================================
fn=flip(conj(hn)); % Matched Filter
fn=fn/sum(fn);
if flag_mod_sim_2==ON
    fn=fn/sqrt(sum(fn.^2));
end
Ef=sum(fn.^2);

%============================================
% Respuesta del Canal Equivalente(rho[n]=h[n]*f[n])
%============================================
rho=conv(hn,fn);

%============================================
% Generacion Simbolos
%============================================
dataIn=randi([0 1],n_symbols,log2(ModOrd));
dataSym=bi2de(dataIn);
ak=transp(qammod(dataSym,ModOrd,0,'gray'));
Ea=mean(abs(ak).^2);
xn = zeros(1,n_symbols*M);
xn(1:M:end) = ak;

%============================================
% Generacion de Componente de Ruido (z[n])
%============================================
No=Ea/SNR;
zn=1.0*randn(1,n_symbols*M)+j*randn(1,n_symbols*M); %ruido complejo Gaussiano
zn=zn/sqrt(var(zn));%potencia 1
zn=sqrt(Eh*No)*zn;

%============================================
% Componente de Senal a la Salida del Canal
%============================================
sn = filter(hn,1,xn);
%break

%============================================
% Senal Rx antes del Filtro Rx
%============================================
rn = sn+zn; % Usar esta senal para el Lab. 4!

%============================================
% Componentes a la Salida del Filtro Rx
%============================================
sn_f = filter(fn,1,sn);
zn_f = filter(fn,1,zn);

%============================================
% Senal Total a la Salida del Filtro Rx 
%============================================
yn = sn_f + zn_f;

%============================================
% Busqueda de Fase Optima
%============================================
min_mse=1e5;
for n=1:M
    sk=sn_f(n:M:20000);
    dataSym_out=qamdemod(sk,ModOrd,0,'gray');
    ak_hat=qammod(dataSym_out,ModOrd,0,'gray');
    error=sk-ak_hat;
    mse=var(error);
    if mse<min_mse
        min_mse=mse;
        phase_op=n;
    end
end
yn=yn(phase_op:end);
sn_f=sn_f(phase_op:end);

%============================================
% Detector
%============================================
yk=yn(1:M:end);
dataSym_out=qamdemod(yk,ModOrd,0,'gray');
ak_hat=qammod(dataSym_out,ModOrd,0,'gray');

%============================================
% Sincronizacion para Calculo de BER
%============================================
td = finddelay(ak_hat,ak);
if td>=0
    ak_hat=[zeros(1,td) ak_hat(1:end)];
    yk=[zeros(1,td) yk(1:end)];
    Rmin=td+1;
    Rmax=length(ak_hat);
else
    ak_hat=[ak_hat(-td+1:end) zeros(1,-td)];
    yk=[yk(-td+1:end) zeros(1,-td)];
    Rmin=1;
    Rmax=length(ak_hat)+td;
end

%============================================
% SNR en el Slicer (uso simbolos Tx)
%============================================
error_slicer=yk-ak(1:length(yk));
error_slicer=error_slicer(1:end-abs(td));
SNR_slicer=Ea/var(error_slicer);

%============================================
% Deteccion y Calculo de BER
%============================================
dataSym_out=qamdemod(ak_hat,ModOrd,0,'gray');
dataOut = de2bi(dataSym_out,log2(ModOrd));
dataIn=dataIn(Rmin:Rmax,:);
dataOut=dataOut(Rmin:Rmax,:);
[numErrors,ber_sim] = biterr(dataIn,dataOut);

%============================================
% Calculo de BER Teorico
%============================================
if ModOrd==2
    ber_teor=qfunc(sqrt(2*SNR));
    ber_teor_sim=qfunc(sqrt(2*SNR_slicer));
elseif ModOrd==4
    ber_teor=qfunc(sqrt(SNR));
    ber_teor_sim=qfunc(sqrt(SNR_slicer));
elseif ModOrd==16
    ber_teor=3/4*qfunc(sqrt(SNR/5));
    ber_teor_sim=3/4*qfunc(sqrt(SNR_slicer/5));
elseif ModOrd==64
    ber_teor=4/6*(1-1/sqrt(ModOrd))*qfunc(sqrt(3/(ModOrd-1)*SNR));
    ber_teor_sim=4/6*(1-1/sqrt(ModOrd))*qfunc(sqrt(3/(ModOrd-1)*SNR_slicer));
else
    'Modulacion no soportada'
    pause
end

%============================================
% Presentacion de Resultados de Simulacion
%============================================
fprintf('\nSNR @ Slicer = %3.2f dB  /  BER Sim.= %1.3e (%1.3e)',10*log10(SNR_slicer),ber_sim,ber_teor_sim);
fprintf('\nSNR Teor  = %3.1f dB  /  BER Teor.= %1.3e\n',SNRdB,ber_teor);

%============================================
% Generacion de Plots
%============================================
eyediagram(real(sn_f(1+M*200:25000)),M,1)
scatterplot(ak);
title('Constelacion Tx')
scatterplot(yk)
%scatterplot(qq(1,:)+j*qq(2,:))
title('Entrada Slicer')
scatterplot(error_slicer)
title('Ruido')
%break
figure
[p1 w1]=psd(rn.*(-1).^[1:length(rn)],1000);
[p2 w2]=psd(sn.*(-1).^[1:length(rn)],1000);
plot(w1-1,10*log10(p1),'b',w2-1,10*log10(p2),'r')
axis([-1 1 -40 10*log10(Ea*Eh)])


%=======================================================================
% Procesamiento de Señales en Tiempo Discreto
% Prof.: Dr. Mario Hueda (mario.hueda@unc.edu.ar)
% 
%Practico Lab. 3
%=======================================================================


clear;
close all;


%============================================
% Generacion de la Respuesta al Impulso
%============================================


fB = 32e9;	% Velocidad de simbolos (baud rate)
T = 1/fB; % Tiempo entre simbolos
M = 8;  %Factor de sobremuestreo
fs = fB*M;	% Sample rate
beta = 0.5001; %Factor de roll-off
EsNodB = 2:10;               % Determina la relación en decibeles entre la energia del símbolo Es y la densidad espectral de potencia de ruido No
EsNo = 10.^(EsNodB/10);
n_symbols = 1000000000;
L = 20;  % 2*L*M+1 es el largo del filtro sobremuestreado
t = [-L:1/M:L]*T;
n_delay_RC_filter = L*M; %Retardo del filtro
rcc = sinc(t/T).*cos(pi*beta*t/T)./(1-4*beta^2*t.^2/T^2);
positive_carrier_frequency = (2.*pi*(fB/2.)*(1+beta)*1j);
negative_carrier_frequency = (-1)*(2.*pi*(fB/2.)*(1+beta)*1j);

SwitchChannel = 'OFF';        % ON: Enciende el canal de ruido gaussiano,
                              % OFF: Apaga el canal de ruido gaussiano

%gn=rcosine(fB,fs,'normal',beta,22); 
% Tx Filter


figure(1)
h = stem(rcc);
title('Respuesta al Impulso del Filtro RC');
%break

%============================================
% Generacion Simbolos
%============================================

for k=1:length(EsNodB)
       
    first_line_symbols = 2*randint(1,n_symbols)-1+j*(2*randint(1,n_symbols)-1);
    second_line_symbols = 2*randint(1,n_symbols)-1+j*(2*randint(1,n_symbols)-1);
    
    upsampled_symbols_line1 = upsample(first_line_symbols,M);  % Agregado de ceros entre los símbolos
    upsampled_symbols_line2 = upsample(second_line_symbols,M);  % Agregado de ceros entre los símbolos
    

    %figure(3)
    %h = plot(ak,'x');
    %xlabel('Real(ak)')
    %ylabel('Imag(ak)');
    %axis([-1.2 1.2 -1.2 1.2])

    %break
    %============================================
    % Señal Banda-Base
    %============================================

    signal_line1 = conv(upsampled_symbols_line1,rcc);
    signal_line2 = conv(upsampled_symbols_line2,rcc);

    figure(2)
    subplot 511
    [pxx1,f] = pwelch(signal_line1,[],[],8192,fs,'centered', 'psd');
    h1 = plot(f,pow2db(pxx1));
    set(h1,'Color','b');
    title('First line PSD');
    subplot 512
    [pxx2,f] = pwelch(signal_line2,[],[],8192,fs,'centered', 'psd');
    h2 = plot(f,pow2db(pxx2));
    set(h2,'Color','r');
    title('Second line psd');

    %============================================
    % Señal Modulada (Analitica)
    %============================================

    n=[1:length(signal_line1)];
    Omega_c1 =(2*pi*(fB/2.)*(1.+beta))/fs;
    Omega_c2 =(-1)*(2*pi*(fB/2.)*(1.+beta))/fs;
     %Portadora arbitraria para generar señal analitica (ojo: depende de M)
    carrier1=exp(j*Omega_c1*n);
    carrier2=exp(j*Omega_c2*n);

    shifted_signal1 = signal_line1.*carrier1;
    shifted_signal2 = signal_line2.*carrier2;

    tx_signal = shifted_signal1+shifted_signal2;


    subplot 513
    [pxx3,f] = pwelch(shifted_signal1,[],[],[],fs,'centered', 'psd');
    s1 = plot(f,pow2db(pxx3));
    xlabel('Frequency (Hz)');
    ylabel('Psd (dB)');
    %legend('pwelch')
    subplot 514
    [pxx4,f] = pwelch(shifted_signal2,[],[],[],fs,'centered', 'psd');
    s2 = plot(f,pow2db(pxx4));
    set(s2,'Color','r');
    xlabel('Frequency (Hz)');
    ylabel('Psd (dB)');
    subplot 515
    [pxx_tx,f] = pwelch(tx_signal,[],[],[],fs,'centered', 'psd');
    tx = plot(f,pow2db(pxx_tx));
    set(tx,'Color','g');
    xlabel('Frequency (Hz)');
    ylabel('Psd (dB)');
    
    %Normalizacion
    vx = var(tx_signal);
    fvx = vx/EsNo(k);
    
    %Ruido AWGN
    if strcmp(SwitchChannel,'ON')
        Noise = randn(1,length(tx_signal));
        rx_signal = tx_signal + (Noise*sqrt(fvx));
    elseif strcmp(SwitchChannel, 'OFF')
        Noise = zeros(1,length(tx_signal));
        rx_signal = tx_signal + Noise;
    else
       display('Opción no valida') 
    end


    rx_baseb_line1 = rx_signal.*carrier2;
    rx_baseb_line2 = rx_signal.*carrier1;

    rx_filtered_line1 = conv(rx_baseb_line1,rcc); 
    rx_filtered_line2 = conv(rx_baseb_line2,rcc); 

    rx_ds_line1_x     = downsample(rx_filtered_line1, M);
    rx_ds_line1       = rx_ds_line1_x(1:length(first_line_symbols));
    rx_ds_line2_x     = downsample(rx_filtered_line2, M);
    rx_ds_line2       = rx_ds_line2_x(1:length(second_line_symbols));


    figure(3)
    subplot 511
    [pxx1,f] = pwelch(rx_signal,[],[],8192,fs,'centered', 'psd');
    h1 = plot(f,pow2db(pxx1));
    set(h1,'Color','g');
    title('Rx PSD');
    subplot 512
    [pxx2,f] = pwelch(rx_baseb_line1,[],[],8192,fs,'centered', 'psd');
    h2 = plot(f,pow2db(pxx2));
    set(h2,'Color','b');
    title('Rx line1 baseband');
    subplot 513
    [pxx3,f] = pwelch(rx_baseb_line2,[],[],[],fs,'centered', 'psd');
    s1 = plot(f,pow2db(pxx3));
    set(s1,'Color','r');
    xlabel('Frequency (Hz)');
    ylabel('Baseband line2 Psd (dB)');
    %legend('pwelch')
    subplot 514
    [pxx4,f] = pwelch(rx_filtered_line1,[],[],[],fs,'centered', 'psd');
    s2 = plot(f,pow2db(pxx4));
    set(s2,'Color','b');
    xlabel('Frequency (Hz)');
    ylabel('Filtered line1 Psd (dB)');
    subplot 515
    [pxx_tx,f] = pwelch(rx_filtered_line2,[],[],[],fs,'centered', 'psd');
    tx = plot(f,pow2db(pxx_tx));
    set(tx,'Color','r');
    xlabel('Frequency (Hz)');
    ylabel('Filtered line 2 Psd (dB)');
    
    
   
    ipHat1 = zeros(1,length(first_line_symbols));
           % demodulation
    l1_re = real(rx_ds_line1); % real
    l1_im = imag(rx_ds_line1); % imaginary
    ipHat1(find(l1_re < 0 & l1_im < 0)) = -1 + -1*j;
    ipHat1(find(l1_re >= 0 & l1_im > 0)) = 1 + 1*j;
    ipHat1(find(l1_re < 0 & l1_im >= 0)) = -1 + 1*j;
    ipHat1(find(l1_re >= 0 & l1_im < 0)) = 1 - 1*j;

    nErr1(k) = size(find([first_line_symbols- ipHat1]),2); % couting the number of errors
        
    ipHat2 = zeros(1,length(first_line_symbols));
           % demodulation
    l2_re = real(rx_ds_line2); % real
    l2_im = imag(rx_ds_line2); % imaginary
    ipHat2(find(l2_re < 0 & l1_im < 0)) = -1 + -1*j;
    ipHat2(find(l2_re >= 0 & l1_im > 0)) = 1 + 1*j;
    ipHat2(find(l2_re < 0 & l1_im >= 0)) = -1 + 1*j;
    ipHat2(find(l2_re >= 0 & l1_im < 0)) = 1 - 1*j;

    nErr2(k) = size(find([second_line_symbols- ipHat2]),2); % couting the number of errors
    
     
end


   


simSer_QPSK = nErr1/length(first_line_symbols);
theorySer_QPSK = erfc(sqrt(0.5*(10.^(EsNodB/10)))) - (1/4)*(erfc(sqrt(0.5*(10.^(EsNodB/10))))).^2;
close all
figure
semilogy(EsNodB,theorySer_QPSK,'b.-');
hold on
semilogy(EsNodB,simSer_QPSK,'mx-');
axis([-3 15 10^-5 1])
grid on
legend('theory-QPSK', 'simulation-QPSK');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for QPSK(4-QAM)')








    %figure(3)
    %subplot 211
    %asd1 = spectrum.welch; 
    %Hpsd_1 = psd(asd1,shifted_signal1,'nfft',8192, 'Fs', fs,'SpectrumType','Onesided');
    %pl1 = plot(Hpsd_1);
    %set(pl1,'Color','b');
    %title('Shifted signal 1');
    %subplot 212
    %asd2 = spectrum.welch; 
    %Hpsd_2 = psd(asd2,shifted_signal2,'nfft',8192, 'Fs', fs, 'SpectrumType','Onesided');
    %pl2 = plot(Hpsd_2)
    %set(pl2,'Color','r');
    %title('Shifted signal 2')






    %break

    %============================================
    % Señal Transmitida (Parte Real de la Señal Analitica)
    %============================================
    %sn_r = real(sn);
    %figure(6)
    %q = spectrum.welch;
    %Hpsd = psd(q,sn_r,'nfft',1024,'SpectrumType','twosided');
    %h=plot(Hpsd);
    %title('PSD Señal Transmitida (Real)')
    %break

    %============================================
    % Filtro de Particion de Fase
    %============================================


    %sn_i = conv(sn_r,fn);
    %sn_hat = sn_r+j*sn_i(Lf+1+0:end-Lf+0);
    %figure(7)
    %q = spectrum.welch;
    %Hpsd = psd(q,sn_hat,'nfft',1024,'SpectrumType','twosided');
    %h=plot(Hpsd);
    %title('PSD Señal a la Salida del Filtro de Particion de Fase')
    %break

    %============================================
    % Comparacion Partes Imaginarias (transmitida y recibida) 
    %============================================
    %sn_i = conv(sn_r,fn);
    %sn_hat = sn_r+j*sn_i(Lf+1:end-Lf);
    %figure(8)
    %n=[1:100]+1000; %100 puntos de ventana de tiempo arbitraria
    %h=plot(n,imag(sn(n)),'b',n,sn_i(n+Lf),'ro');
    %legend('Transmitida', 'Recuperada');
    %title('Parte Imaginaria de la Señal Analitica')
    %break
    %============================================
    % Señal Demodulada (Banda-Base)
    %============================================

    %n=[1:length(sn_hat)];
    %carrier=exp(j*Omega_c*n);
    %yn_hat = sn_hat.*conj(carrier); 
    %portadora conjugada (demodulacion)
    %figure(9)
    %q = spectrum.welch;
    %Hpsd = psd(q,yn_hat,'nfft',1024,'SpectrumType','twosided');
    %h=plot(Hpsd);
    %title('PSD Señal Demodulada (Banda-Base)')


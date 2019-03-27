%=======================================================================
% Comunicaciones Opticas Coherentes
% Prof.: Dr. Mario Hueda (mario.hueda@unc.edu.ar)
% Practico Lab. 1
%=======================================================================
clear;
close all;

%============================================
% Generacion de la Respuesta al Impulso
%============================================
fB = 32e9;	% Velocidad de simbolos (baud rate)
T = 1/fB; % Tiempo entre simbolos
M = 16;  %Factor de sobremuestreo
fs = fB*M;	% Sample rate
ModOrd=2^4;
SNRdB=[1:10];
SNR=10.^(SNRdB/10);
n_symbols = 500000;
flag_mod_sim_2=0; % "0" Modelo 1 / "1" Modelo 2

%============================================
% Respuesta del Filtro Transmisor
%============================================
beta = 0.25; %Factor de roll-off
gn=rcosine(fB,fs,'normal',beta,24);  % Generacion usando funcion de matlab 
gn=rcosine(fB,fs,'sqrt',beta,24);  %  Generacion usando funcion de matlab (raiz cuadrada)
BW_EF=fB*.4;
[he_b,he_a] = butter(3,2*BW_EF/(M*fB));
%gn= filter(he_b',he_a',[0 0 0 1 zeros(1,200)]);

%Obtencion de energia
gn=gn/sum(gn)*M;
Eg=sum(gn.^2);
%break

figure(1)
h = stem(gn);
title('Respuesta al Impulso g[n]');
xlabel('n');
grid

%============================================
% Respuesta del Filtro Transmisor+Canal (h[n])
%============================================
bn=1;
hn=conv(gn,bn);
hn=hn/sum(hn)*M;

if flag_mod_sim_2==1
    hn=hn/sqrt(sum(hn.^2));
end
Eh=sum(hn.^2);
%break

figure(2)
h = stem(hn);
title('Respuesta al Impulso h[n]');
xlabel('n');
grid

%break

%============================================
% Respuesta del Filtro Receptor (f[n])
%============================================
fn=flip(conj(hn)); % Matched Filter
BW_EF=fB*.5;
[he_b,he_a] = butter(3,2*BW_EF/(M*fB));
%fn= filter(he_b',he_a',[0 0 0 1 zeros(1,100)]);

fn=fn/sum(fn);
if flag_mod_sim_2==1
    fn=fn/sqrt(sum(fn.^2));
end
Ef=sum(fn.^2);


figure(3)
h = stem(fn);
title('Respuesta al Impulso f[n]');
xlabel('n');
grid

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
%break

figure(4)
h = stem(xn(10:10+M*10));
ylabel('x[n]')
xlabel('n');

for k=1:length(SNRdB)
        %============================================
        % Generacion de Componente de Ruido (z[n])
        %============================================
        No=Ea/SNR(k);
        zn=1/sqrt(2)*(randn(1,n_symbols*M)+j*randn(1,n_symbols*M)); %potencia 1
        zn=sqrt(Eh*No)*zn;
        %break

        %============================================
        % Componente de Senal a la Salida del Canal
        %============================================
        sn = filter(hn,1,xn);
        %break

        eg=0.1;
        phi_e=1.0*pi/8;
        k1=0.5*((1-eg)*exp(j*phi_e/2)+(1+eg)*exp(-j*phi_e/2));
        k2=0.5*((1-eg)*exp(j*phi_e/2)-(1+eg)*exp(-j*phi_e/2));
        sn= k1*sn+k2*conj(sn);
        k1I=real(k1);k2I=real(k2);
        k1Q=imag(k1);k2Q=imag(k2);
        W=[k1I+k2I  -k1Q+k2Q;k1Q+k2Q k1I-k2I];
        Winv=inv(W);
        %Winv=eye(2,2);
        akd=W(1,1)*real(ak)+W(1,2)*imag(ak)+j*(W(2,1)*real(ak)+W(2,2)*imag(ak));
        %break

        %var(sn)
        %break
        %============================================
        % Componentes a la Salida del Filtro Rx
        %============================================
        sn_f = filter(fn,1,sn);
        zn_f = filter(fn,1,zn);
        %break

        %============================================
        % Senal Total a la Salida del Filtro Rx 
        %============================================
        yn = sn_f + zn_f;

        %============================================
        % SNR
        %============================================

        SNR_sim(k)=Ea*max(rho)^2/(Eh*Ef*No);

       
        min_mse=1e5;
        for n=1:M
            sk=sn_f(n:M:20000);
            yk=yn(n:M:20000);
            sk=Winv(1,1)*real(sk)+Winv(1,2)*imag(sk)+j*(Winv(2,1)*real(sk)+Winv(2,2)*imag(sk));
            yk=Winv(1,1)*real(yk)+Winv(1,2)*imag(yk)+j*(Winv(2,1)*real(yk)+Winv(2,2)*imag(yk));
            ak_hat=my_qamdemod(sk,ModOrd);
            error=yk-ak_hat;
            mse=var(error);
            if mse<min_mse
                min_mse=mse;
                phase_op=n;
            end
        end
        yn=yn(phase_op:end);
        sn_f=sn_f(phase_op:end);

        sn_f=Winv(1,1)*real(sn_f)+Winv(1,2)*imag(sn_f)+j*(Winv(2,1)*real(sn_f)+Winv(2,2)*imag(sn_f));
        yn=Winv(1,1)*real(yn)+Winv(1,2)*imag(yn)+j*(Winv(2,1)*real(yn)+Winv(2,2)*imag(yn));

        %============================================
        % Generacion de Diagrama Ojo
        %============================================
        eyediagram(real(sn_f(1+M*50:5000)),M,1)

        %============================================
        % Detector
        %============================================
        yk=yn(1:M:end);
        ak_hat=my_qamdemod(yk,ModOrd);

        %============================================
        % Sincronizacion y Calculo de BER
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
        error=yk-ak(1:length(yk));
        SNR_slicer=Ea/var(error);
        10*log10(SNR_slicer)
        %break

        dataSym_out=qamdemod(ak_hat,ModOrd,0,'gray');
        dataOut = de2bi(dataSym_out,log2(ModOrd));

        dataIn=dataIn(Rmin:Rmax,:);
        dataOut=dataOut(Rmin:Rmax,:);
        [numErrors,ber_sim] = biterr(dataIn,dataOut);

    %============================================
    % Calculo de BER Teorico
    %============================================
    if ModOrd==2
        ber_teor(k)=qfunc(sqrt(2*SNR(k)));
        ber_teor(k)=qfunc(sqrt(2*SNR_sim(k)));
    elseif ModOrd==4
        ber_teor(k)=qfunc(sqrt(SNR(k)));
        ber_teor_sim(k)=qfunc(sqrt(SNR(k)));
    elseif ModOrd==16
        ber_teor(k)=3/4*qfunc(sqrt(SNR(k)/5));
        ber_teor(k)=3/4*qfunc(sqrt(SNR_sim(k)/5));
    elseif ModOrd==64
        ber_teor(k)=4/6*(1-1/sqrt(ModOrd))*qfunc(sqrt(3/(ModOrd-1)*SNR(k)));
        ber_teor_sim(k)=4/6*(1-1/sqrt(ModOrd))*qfunc(sqrt(3/(ModOrd-1)*SNR_sim(k)));
    else
        'Modulacion no soportada'
        pause
    end

   
end

figure
semilogy(SNRdB,ber_teor,'b.-');
hold on
semilogy(SNRdB, ber_sim,'mx-');
grid on
legend('theory-QPSK', 'simulation-QPSK');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for QPSK(4-QAM)')






%[ber_sim ber_teor_sim ber_teor]

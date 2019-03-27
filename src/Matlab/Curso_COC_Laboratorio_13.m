%=======================================================================% Comunicaciones Opticas Coherentes% Prof.: Dr. Mario Hueda (mario.hueda@unc.edu.ar)% Practico Laboratorio%=======================================================================clear;close allN=74321234;rand('seed',N);randn('seed',round(sqrt(N)));ON=1;OFF=0;%-----------------------------------------------------------------------% Parametros Generales del Sistema%-----------------------------------------------------------------------M=4;                % Oversampling Factorf_b = 32*1;			% Symbol data rate [GBd]f_s = f_b*M;   	% Sampling rate [GHz]n_data=1000000;n_data_cma=round(0.2*n_data); % Poner en 0 cuando solo se use LMSModOrd=2^4;OSNRdB=20; %NOTA: En Laboratorio 16 aumentar 3dB para 64GBdBW_EX_Tx=10;     %BW Excess [%]BW_EX_Rx=BW_EX_Tx; Dv_RX =.50*.0;  %Linewidth in MHzcarrier_offset=.00; %GHz (En Laboratorio 16, el rango de captura es ~ 100MHz)%-----------------------------------------------------------------------% Calculo de SNR%-----------------------------------------------------------------------Optical_BW=f_b; %Baud Rate per PolarizationOptical_BWref=12.5;%0.1nmSNR_correction_factor_dB=10*log10(Optical_BWref/Optical_BW);SNRdB=OSNRdB+SNR_correction_factor_dB;%-----------------------------------------------------------------------% Parametros del Canal Opticos%-----------------------------------------------------------------------DGD=0;    %[ps]SOPMD=000; %[ps2]SOPMD_PCD=1000; %[ps2]L=0;       %[km]beta2=21.67; %ps2Leq=SOPMD_PCD/2/beta2; %[km]                       %L=-Leq;       %[km]%-----------------------------------------------------------------------% Show Simulation Setup%-----------------------------------------------------------------------fprintf('\nBaud Rate=%3.1f GBd - BW Excess=%2.0f - OSNR=%2.2f dB (%2.2f dB)',f_b,BW_EX_Tx,OSNRdB,SNRdB);%=======================================================================%   Generacion de Simbolos%=======================================================================n_data_bits=n_data;dataIn_x=randi([0 1],n_data_bits,log2(ModOrd));dataSym_x=bi2de(dataIn_x);symb_txd_x=transp(qammod(dataSym_x,ModOrd,0,'gray'))*1.;sb_x= zeros(1,f_s*n_data/f_b);sb_x(1:f_s/f_b:end)=symb_txd_x;dataIn_y=randi([0 1],n_data_bits,log2(ModOrd));dataSym_y=bi2de(dataIn_y);symb_txd_y=transp(qammod(dataSym_y,ModOrd,0,'gray'));sb_y= zeros(1,f_s*n_data/f_b);sb_y(1:f_s/f_b:end)=symb_txd_y;Ea=var(symb_txd_x); %Potencia de la constelacion por polarizacion%=======================================================================%	Filtro Tx%=======================================================================htxd=rcosine(f_b,f_s,'sqrt',BW_EX_Tx/100,22);  % Tx Filterhtxd=htxd/sqrt(sum(htxd.^2));%=======================================================================%	Filtro Rx%=======================================================================hrxd=rcosine(f_b,f_s,'sqrt',BW_EX_Rx/100,22);  % Tx Filterhrxd=hrxd/sqrt(sum(hrxd.^2));%=======================================================================%	Generacion de la Respuesta del Canal MIMO%=======================================================================[h11 h12 h21 h22]=generador_fibra_cd_pmd(DGD,SOPMD,SOPMD_PCD,L,M,f_b,htxd,1);%=======================================================================% Senal Optica Recibida%=======================================================================s_x=filter(h11,1,sb_x)+filter(h12,1,sb_y);s_y=filter(h21,1,sb_x)+filter(h22,1,sb_y);%=======================================================================% Laser Phase Noise (@ Rx)%=======================================================================t_line = [1/f_s:1/f_s:n_data/f_b]; % Vector con muestras de tiempo a 1/f_sDv_RX = Dv_RX/1e3;sigma_lo = sqrt(2*pi*Dv_RX/f_s);np=sigma_lo*randn(1,length(s_x));t_rxd=t_line(1:length(s_x));laser_phase_noise=cumsum(np);%=======================================================================% Generacion Ruido ASE%=======================================================================No=10^(-SNRdB/10)*Ea;n_x=randn(1,length(s_x))+j*randn(1,length(s_x));n_y=randn(1,length(s_y))+j*randn(1,length(s_y));n_x=sqrt(No/2)*n_x; %Potencia Non_y=sqrt(No/2)*n_y; %Potencia No%=======================================================================% Senal Optica Recibida%=======================================================================p_x=s_x+n_x;p_y=s_y+n_y;p_x=p_x.*exp(j*2*pi*carrier_offset*t_rxd+j*laser_phase_noise);p_y=p_y.*exp(j*2*pi*carrier_offset*t_rxd+j*laser_phase_noise);%=======================================================================% Filtrado Rx%=======================================================================r_x=filter(hrxd,1,p_x);r_y=filter(hrxd,1,p_y);%=======================================================================% Simulador del Receptor%=======================================================================OS=2; %Sobremuestreo del RxP=16*1.; % Factor de paralelismoP_OS=round(P*OS);kp=0.5; % Ganancia proporcional PLL ki=kp/200; % Ganancia integral PLLNFFE=16*1; % Nro de taps del ecualizador T/2t_rxd_2baud=[1/f_s:1/OS/f_b:length(r_x)/f_s];srxd_2baud_x=interp1(t_line,r_x,t_rxd_2baud,'linear');srxd_2baud_y=interp1(t_line,r_y,t_rxd_2baud,'linear');N1=NFFE/2;Coef_FFE=zeros(1,NFFE);Coef_FFE(N1)= 1.;Coef_FFE_xx=Coef_FFE;Coef_FFE_xy=Coef_FFE*0;Coef_FFE_yx=Coef_FFE*0;Coef_FFE_yy=Coef_FFE;nco_x=0;nco_y=0;ko=250;symb_rxd_x=zeros(1,n_data+ko);symb_rxd_y=zeros(1,n_data+ko);index_k=NFFE;sum_phe_x=0;sum_phe_y=0;ffe_dd_start=OFF;fcr_start=OFF;mse_x=0;mse_y=0;show_plot=400;show_plot_index=1;phi0=mean(abs(srxd_2baud_x).^2);R2=(mean(abs(symb_txd_x).^4)/mean(abs(symb_txd_x).^2)); %Potencia de la constelacion por polarizacionclose allfor k=1:P:n_data-250    %==================================================    %    FSE (OS)    %==================================================    Input_Signal_x=srxd_2baud_x(index_k:index_k+NFFE+P_OS-2);    Input_Signal_y=srxd_2baud_y(index_k:index_k+NFFE+P_OS-2);    m_xx=convmtx(Coef_FFE_xx,P_OS);    m_xy=convmtx(Coef_FFE_xy,P_OS);    m_yx=convmtx(Coef_FFE_yx,P_OS);    m_yy=convmtx(Coef_FFE_yy,P_OS);    Output_x_To2=m_xx*transp(Input_Signal_x)+m_xy*transp(Input_Signal_y);    Output_y_To2=m_yx*transp(Input_Signal_x)+m_yy*transp(Input_Signal_y);    %==================================================    %    Baud Rate Convertion (OS/T-->1/T) and     %    Carrier Recovery    %==================================================    Output_x_baud=Output_x_To2(1:2:end)*exp(-j*nco_x);    Output_y_baud=Output_y_To2(1:2:end)*exp(-j*nco_y);        %==================================================    %    Slicer    %==================================================    out_slicer_x=transp(my_qamdemod(Output_x_baud,ModOrd));    out_slicer_y=transp(my_qamdemod(Output_y_baud,ModOrd));    symb_rxd_x(k+ko:k+ko+P-1)=out_slicer_x;        symb_rxd_y(k+ko:k+ko+P-1)=out_slicer_y;        %==================================================    %    State Machine    %==================================================    if abs(k-(n_data_cma+0000))<P        ffe_dd_start=ON;        fcr_start=ON;    end        %==================================================    %    Error Computation and Interpolation    %==================================================    if (ffe_dd_start)        error_x=transp(out_slicer_x-transp(Output_x_baud));        error_y=transp(out_slicer_y-transp(Output_y_baud));    else        error_x=Output_x_baud.*(sqrt(R2)-abs(Output_x_baud));        error_y=Output_y_baud.*(sqrt(R2)-abs(Output_y_baud));    end    error_xx=transp(transp(Output_x_baud)-out_slicer_x);    error_yy=transp(transp(Output_y_baud)-out_slicer_y);    mse_x=.999*mse_x+.001*mean(error_xx.*conj(error_xx));    mse_y=.999*mse_y+.001*mean(error_yy.*conj(error_yy));    % Modulacion del error    error_x=error_x*exp(j*nco_x);    error_y=error_y*exp(j*nco_y);        %==================================================    %    LMS Algorithm    %==================================================    if (ffe_dd_start)        step_fff=1*.25/(1*NFFE*phi0*P);    else        step_fff=.25/(5*NFFE*phi0*P);    end    accum_error_signal_xx=zeros(1,NFFE);    accum_error_signal_xy=zeros(1,NFFE);    accum_error_signal_yx=zeros(1,NFFE);    accum_error_signal_yy=zeros(1,NFFE);    for l=1:1:P        accum_error_signal_xx=accum_error_signal_xx+error_x(l)*conj(Input_Signal_x(2*l-1:NFFE+2*l-2));        accum_error_signal_xy=accum_error_signal_xy+error_x(l)*conj(Input_Signal_y(2*l-1:NFFE+2*l-2))*1.0;        accum_error_signal_yx=accum_error_signal_yx+error_y(l)*conj(Input_Signal_x(2*l-1:NFFE+2*l-2))*1.0;        accum_error_signal_yy=accum_error_signal_yy+error_y(l)*conj(Input_Signal_y(2*l-1:NFFE+2*l-2));    end    %----------------------------------------------    % Taps update    %----------------------------------------------    Coef_FFE_xx=Coef_FFE_xx+step_fff*accum_error_signal_xx;    Coef_FFE_xy=Coef_FFE_xy+step_fff*accum_error_signal_xy;    Coef_FFE_yx=Coef_FFE_yx+step_fff*accum_error_signal_yx;    Coef_FFE_yy=Coef_FFE_yy+step_fff*accum_error_signal_yy;        %----------------------------------------------    % Singularity Prevent Procedure (default: turned off)    %----------------------------------------------    if abs(k-(n_data_cma/2))<-P         Coef_FFE_yy=conj(Coef_FFE_xx(end:-1:1));        Coef_FFE_yx=-conj(Coef_FFE_xy(end:-1:1));        display('Pol. Colapse Prevent Done');    end        %==================================================    %    Fine Carrier Recovery    %==================================================    if (fcr_start)        phe_x=mean(angle(conj(out_slicer_x).*transp(Output_x_baud)));        phe_y=mean(angle(conj(out_slicer_y).*transp(Output_y_baud)));        sum_phe_x=phe_x+sum_phe_x;        sum_phe_y=phe_y+sum_phe_y;        nco_x=nco_x+kp*phe_x+ki*sum_phe_x;        nco_y=nco_y+kp*phe_y+ki*sum_phe_y;    end            %==================================================    %   Plots     %==================================================    srxd_baudn_x(k:k+P-1)=transp(Output_x_baud);    srxd_baudn_y(k:k+P-1)=transp(Output_y_baud);    index_k=index_k+P_OS;        if show_plot_index==show_plot*1        show_plot_index=1;        subplot 121        plot([1:NFFE],real(Coef_FFE_xx),'rx-',[1:NFFE],imag(Coef_FFE_xx),'ro-');        hold on        plot([1:NFFE],real(Coef_FFE_xy),'rx:',[1:NFFE],imag(Coef_FFE_xy),'ro:');        plot([1:NFFE],real(Coef_FFE_yx),'bx:',[1:NFFE],imag(Coef_FFE_yx),'bo:');        plot([1:NFFE],real(Coef_FFE_yy),'bx-',[1:NFFE],imag(Coef_FFE_yy),'bo-');        title(sprintf('Percentage = %3.0f',round(k/(n_data-250)*100)));        hold off                subplot 122        if k>1000            plot(srxd_baudn_x(k-1000:k),'rx');hold on;plot(srxd_baudn_y(k-1000:k),'bo');hold off            axis(1.5*[-log2(ModOrd) log2(ModOrd) -log2(ModOrd) log2(ModOrd)])            grid           hold off           title(sprintf('SNR = %3.2f/%3.2f dB',10*log10(Ea/mse_x),10*log10(Ea/mse_y)));        end        drawnow    end    show_plot_index=show_plot_index+1;endxc_xx=xcorr(symb_txd_x(end-5000:end),symb_rxd_x(end-5000:end));xc_yx=xcorr(symb_txd_y(end-5000:end),symb_rxd_x(end-5000:end));if max(abs(xc_xx))>max(abs(xc_yx))    [ser_x,ber_x]=ber_ser_comp(symb_txd_x(1:end-100),symb_rxd_x(1:end-100-ko),n_data_cma+150000,ModOrd);    [ser_y,ber_y]=ber_ser_comp(symb_txd_y(1:end-100),symb_rxd_y(1:end-100-ko),n_data_cma+150000,ModOrd);else    [ser_x,ber_x]=ber_ser_comp(symb_txd_x(1:end-100),symb_rxd_y(1:end-100-ko),n_data_cma+150000,ModOrd);    [ser_y,ber_y]=ber_ser_comp(symb_txd_y(1:end-100),symb_rxd_x(1:end-100-ko),n_data_cma+150000,ModOrd);end    fprintf('\nSER=%f / %f',ser_x, ser_y);fprintf('\nBER=%f / %f',ber_x, ber_y);fprintf('\nSNRx/SNRy (Averg.)=%3.2f / %3.2f  (%3.2f)',10*log10(Ea/mse_x),10*log10(Ea/mse_y),5*log10(Ea/mse_x)+5*log10(Ea/mse_y));fprintf('\n');[Q(sqrt(10^(SNRdB/10))) 3/4*Q(sqrt(10^(SNRdB/10)/5))]
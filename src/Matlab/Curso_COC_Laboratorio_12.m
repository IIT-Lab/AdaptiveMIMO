%=======================================================================% Comunicaciones Opticas Coherentes% Prof.: Dr. Mario Hueda (mario.hueda@unc.edu.ar)% Practico Laboratorio%=======================================================================clear;close allN=74321234;rand('seed',N);randn('seed',round(sqrt(N)));M=8;                % Oversampling Factorf_b = 32*1.;			% Symbol data ratef_s = f_b*M;   	% Sampling raten_data=100000;ModOrd=2^4;BW_EX_Tx=10;     %BW Excess [%]BW_EX_Rx=BW_EX_Tx;     %BW Excess [%]%-----------------------------------------------------------------------% Parametros del Canal Optico%-----------------------------------------------------------------------DGD=100;    %[ps]SOPMD=4000; %[ps2]SOPMD_PCD=750; %[ps2]L=500; %[km]%=======================================================================%   Generacion de Simbolos%=======================================================================n_data_bits=n_data;dataIn_x=randi([0 1],n_data_bits,log2(ModOrd));dataSym_x=bi2de(dataIn_x);symb_txd_x=transp(qammod(dataSym_x,ModOrd,0,'gray'))*1.;sb_x= zeros(1,f_s*n_data/f_b);sb_x(1:f_s/f_b:end)=symb_txd_x;dataIn_y=randi([0 1],n_data_bits,log2(ModOrd));dataSym_y=bi2de(dataIn_y);symb_txd_y=transp(qammod(dataSym_y,ModOrd,0,'gray'));sb_y= zeros(1,f_s*n_data/f_b);sb_y(1:f_s/f_b:end)=symb_txd_y;Ea=var(symb_txd_x); %Potencia de la constelacion%=======================================================================%	Filtro Tx%=======================================================================htxd=rcosine(f_b,f_s,'sqrt',BW_EX_Tx/100,22);  % Tx Filterhtxd=htxd/sqrt(sum(htxd.^2));%=======================================================================%	Filtro Rx%=======================================================================hrxd=rcosine(f_b,f_s,'sqrt',BW_EX_Rx/100,22);  % Tx Filterhrxd=hrxd/sqrt(sum(hrxd.^2));%=======================================================================%	Generacion de la Respuesta del Canal MIMO%=======================================================================[h11 h12 h21 h22]=generador_fibra_cd_pmd(DGD,SOPMD,SOPMD_PCD,L,M,f_b,htxd,1);%=======================================================================% Mostrar Setup de Simulacion%=======================================================================fprintf('Baud Rate=%3.1f GBd - BW Excess=%2.0f - L=%3.0f km- DGD=%3.0f ps- SOPMD=%4.0f ps2\n',f_b,BW_EX_Tx,L, DGD,SOPMD);%break%=======================================================================% Senal Optica Recibida%=======================================================================s_x=filter(h11,1,sb_x)+filter(h12,1,sb_y);s_y=filter(h21,1,sb_x)+filter(h22,1,sb_y);%=======================================================================% Filtro Apareado%=======================================================================q_x=filter(conj(flip(h11)),1,s_x)+filter(conj(flip(h21)),1,s_y);q_y=filter(conj(flip(h12)),1,s_x)+filter(conj(flip(h22)),1,s_y);%=======================================================================% Generacion Diagrama Ojos%=======================================================================eyediagram(s_x(M*(100):M*(200)),M,1,M/2)eyediagram(q_x(M*(100):M*(200)),M,1,M/2)
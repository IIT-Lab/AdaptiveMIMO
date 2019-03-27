function [h11 h12 h21 h22]=generador_fibra_cd_pmd(DGD,SOPMD,SOPMD_PCD,L,R,Br,htxd,PLOT)
% DGD in ps
% SOPMD and SOPMD_PCD in ps^2
% Br in GHz
% L in km

% Fiber parameters:
D=17*1.0;                                       % Dispersion parameter, in ps/(nm*km)
S=0.09*1.0;                                     % Slope parameter, in ps/(km*nm^2)
lambda=1550;                                % In nm
beta2=D*(1e-11)*(lambda^2)/(6*pi);          % In units of 1/(Gigaradian^2 * km / s^2)
beta3=(S*1e-21*lambda^4)/(36*pi^2)- ...
    (lambda*beta2*1e-8)/(3*pi);             % In units of 1/(Gigaradian^3 * km / s^3)

% Other transmission parameters:
fB=Br;                                    % Symbol rate, in GHz
fs=fB*R;                                    % In GHz
fftsize=2^15;                             % Number of points in FFT
fftsize2=fftsize/2;
fstep=fs/fftsize;                           % Frequency bins in the FFT
fgrid=fstep*[0:fftsize2-1,-fftsize2:-1];    % Negative frequencies are shifted to the region [fftsize2:fftsize-1]
omegagrid=2*pi*fgrid;                       % The frequency mask in gigaradians/sec
W=omegagrid;


% Create chromatic dispersion filtering mask 
% to be applied to data as a filter in the frequency domain:
Dl=i*beta2*(omegagrid.^2)/2+i*beta3*(omegagrid.^3)/6;

% Filter pulse in frequency domain:
fft_A=fft(htxd,fftsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          BRUYERE MODEL                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SOPMD=SOPMD*1e-6;
DGD=DGD*1e-3; %DGD in ns
tw=SOPMD_PCD*1e-6; %PCD component in ns^2 

% Si la SOPMD total es menor que PCD--> PCD=Despolarizacion
if SOPMD<tw
    tw=.707*SOPMD;
end
pw_DGD=sqrt(SOPMD.^2-tw^2);

if DGD>0
    kw=.25*pw_DGD/DGD;
else
    kw=0;
end

kw_W=kw*W*1.0;
R1=[cos(kw_W)  -sin(kw_W); 
    sin(kw_W)   cos(kw_W)];
R1t=[cos(kw_W)  sin(kw_W); 
   -sin(kw_W)   cos(kw_W)];
D=[exp(-j*(DGD+tw*W/2).*W/2)  zeros(1,fftsize); 
   zeros(1,fftsize)  exp(j*(DGD+tw*W/2).*W/2)];
DR1t=product_freq_matrices(D,R1t,fftsize);
U=product_freq_matrices(R1,DR1t,fftsize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          DGD AND SOPMD CALCULATION SECTION                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get components of Jones matrix:
U11=U(1,1:fftsize);
U12=U(1,fftsize+1:end);
U21=U(2,1:fftsize);
U22=U(2,fftsize+1:end);


% Compute derivative of Jones matrix:
ks=1.0/(2.0*pi*fstep);
dU11 = ks*([U11(2:end),U11(1)]-U11);
dU12 = ks*([U12(2:end),U12(1)]-U12);
dU21 = ks*([U21(2:end),U21(1)]-U21);
dU22 = ks*([U22(2:end),U22(1)]-U22);

% Compute product -2*i*(dU/domega)*herm(U):
PO=-2*i*[dU11.*conj(U11)+dU12.*conj(U12), dU11.*conj(U21)+dU12.*conj(U22);...
    dU21.*conj(U11)+dU22.*conj(U12), dU21.*conj(U21)+dU22.*conj(U22)];

% Compute vector Omega (polarization vector in  Stokes space) at the central frequency:
Omega=[real(PO(2,1)),imag(PO(2,1)),real(PO(1,1))];

%[norm(Omega) norm([real(PO(2,121)) imag(PO(2,121)) real(PO(1,121))])]

% Compute derivative of Omega with respect to angular frequency at the central frequency:
dOmega=ks*[real(PO(2,2))-real(PO(2,1)),imag(PO(2,2))-imag(PO(2,1)),real(PO(1,2))-real(PO(1,1))];

% Compute DGD and SOPMD at the central frequency:
DGDt=1.0e3*norm(Omega);               % In picoseconds
SOPMDt=1.0e6*norm(dOmega);            % In squared picoseconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CHANNEL RESPONSE CALCULATION SECTION                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uc=U.*[fft_A.*exp(L*Dl) fft_A.*exp(L*Dl); fft_A.*exp(L*Dl) fft_A.*exp(L*Dl)];

% Find time domain representation of channel matrix:
Uc11=Uc(1,1:fftsize);
Uc12=Uc(1,fftsize+1:end);
Uc21=Uc(2,1:fftsize);
Uc22=Uc(2,fftsize+1:end);
u11=ifft(Uc11,fftsize);
u12=ifft(Uc12,fftsize);
u21=ifft(Uc21,fftsize);
u22=ifft(Uc22,fftsize);

u11=ifftshift(ifft(Uc11,fftsize));
u12=ifftshift(ifft(Uc12,fftsize));
u21=ifftshift(ifft(Uc21,fftsize));
u22=ifftshift(ifft(Uc22,fftsize));


% Get components:
hhr=real(u11);
hhi=imag(u11);
vvr=real(u22);
vvi=imag(u22);
hvr=real(u12);
hvi=imag(u12);
vhr=real(u21);
vhi=imag(u21);

% Para acortar el largo de los pulsos guarda el 99% de la energia
total=abs(u11)+abs(u12)+abs(u21)+abs(u22);
acumulado = cumsum(total.^2);
total = acumulado(end);
usar = find((acumulado>total.*0.0001)&(acumulado<total.*0.9999));
pulsos=[hhr' hhi' hvr' hvi' vhr' vhi' vvr' vvi'];
smfpulse=pulsos(usar,:);
[h11 h12 h21 h22]=imp_resp(smfpulse,R);
if PLOT>0
    fprintf('Response Length=%5.0f [Bauds]',length(usar)/R);    
    disp(sprintf('\nEstimated DGD = %3.2f ps Estimated SOPMD = %5.2f ps^2',DGDt,SOPMDt));    t=[1:length(h11)]/R;
    subplot(2,2,1)
    plot(t,real(h11),'r',t,imag(h11),'b');
    subplot(2,2,2)
    plot(t,real(h12),'r',t,imag(h12),'b');
    subplot(2,2,3)
    plot(t,real(h21),'r',t,imag(h21),'b');
    subplot(2,2,4)
    plot(t,real(h22),'r',t,imag(h22),'b');
    drawnow
end

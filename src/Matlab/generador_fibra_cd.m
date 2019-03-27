function [h11 h12 h21 h22]=generador_fibra_cd(L,R,Br,htxd,PLOT)
% L in km
% Br in GHz


% Fiber parameters:
D=17*1.0;                                       % Dispersion parameter, in ps/(nm*km)
S=0.09*1.0;                                     % Slope parameter, in ps/(km*nm^2)
lambda=1550;                                % In nm
beta2=D*(1e-11)*(lambda^2)/(6*pi);          % In units of 1/(Gigaradian^2 * km / s^2)
beta3=1.0*(S*1e-21*lambda^4)/(36*pi^2)- ...
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

Uc=[fft_A.*exp(L*Dl) fft_A*0; fft_A*0 fft_A.*exp(L*Dl)]; %Solo CD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          CHANNEL RESPONSE CALCULATION SECTION                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    fprintf('Response Length=%5.0f [Bauds]\n',length(usar)/R);    
    t=[1:length(h11)]/R;
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

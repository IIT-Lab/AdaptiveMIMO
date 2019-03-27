function [ser,ber]=ber_ser_comp(Atx,Ahat,start,M);
N=length(Atx);
xc=xcorr(Atx,Ahat);
m=find(abs(xc)==max(abs(xc)));
if N>m
    Atx=[zeros(1,N-m) Atx(1:end-N+m)];
else
    Atx=[Atx(m-N+1:end) zeros(1,m-N)];
end
p=round(angle(mean((Atx.*conj(Ahat))))/(pi/2));
Ahat_rot=Ahat*exp(j*p*pi/2);
aux=Atx.*conj(round(Ahat_rot));
m=find(abs(angle(aux(start:end-1000)))>0.05);
ser=length(m)/(N-start-1000);
aux=Atx-round(Ahat_rot);
m=find(aux(start:end-10000)~=0);
ser=length(m)/(N-start-10000);



txSym = qamdemod(Atx(start:end-10000),M,0,'gray');
dataIn = de2bi(txSym,log2(M));
rxSym = qamdemod(round(Ahat_rot(start:end-10000)),M,0,'gray');
dataOut = de2bi(rxSym,log2(M));
nErrors = biterr(dataIn,dataOut);
ber=nErrors/length(dataIn)/log2(M);
function [s, z] = my_qamdemod(y, m)

%this is the function for demodulation of QAM, y is the M-QAM Modulated
%signam, m is playing the same role. z is the demodulated output. input y
%is a row or coulmn vector and, z should be the row vector.

%santosh shah, The LNM IIT Jaipur (India)(santosh.jnt@gmail.com) 23/04/07

% again i need to check the value of M, so by the same process.
if log2(sqrt(m))~= floor(log2(sqrt(m)))
    %error('Please check the value of m that you have provided for type M-QAM.');
end

%taking the reverse process
k = sqrt(m);
r = 2*(0:k-1) - k + 1;
[xi, yi] = meshgrid(r);
c = xi + j*flipud(yi);
c = c(:);

%now comparing the data from c's vector after rounding the input data.
% Allocate space for output
z = zeros(size(y));
s = zeros(size(y));

% Slicer: Find closest constellation symbol, symbol-by-symbol.
for k = 1:length(y)
    [nil ind] = min(abs(y(k) - c));
    z(k) = ind - 1;
    s(k) = c(ind);
end

if m==2
    s=ones(size(y));
    n=find(y<0);
    s(n)=-1;
    %s=s*exp(-j*pi/4);
end
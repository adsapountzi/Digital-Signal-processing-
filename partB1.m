%Author: Sapountzi Athanasia Despoina 02624
%%MEROS B1

%% B1.1
%signal x[n] = 80δ[n] - 80sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16)
%delta -> dirac function

N = 10000;
n= 0:N-1;
Nbig = 10000;
Nbins = 8000;
signal = 80*dirac(n) - 80*sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16);
L =8;
%rectangular window
wrect = rectwin(L);
%apply windowing technique to signal
for n = 1:8
    yrect(n) = (80*dirac(n) - 80*sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16))*wrect(n);
end
figure(1);
plot(yrect);
%% B1.2
N = L+1;
Xk = fftshift(fft(yrect,N));
figure(3);
plot(20*log10(abs(Xk)));

N = 10000;
Xk = fftshift(fft(yrect,N));
figure(4);
plot(20*log10(abs(Xk)));
%% B1.3
L =16;
%rectangular window
wrect = rectwin(L);
%apply window to signal
for n = 1:16
    yrect(n) = (80*dirac(n) - 80*sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16))*wrect(n);
end
%plot( 20*log10(abs(fftshift(fft(yrect, Nbins)))), 'm' ); 
N = L+1;
Xk = fftshift(fft(yrect,N));
figure(6);
plot(20*log10(abs(Xk)));

N = 10000;
Xk = fftshift(fft(yrect,N));
figure(7);
plot(20*log10(abs(Xk)));

L = 32;
%rectangular window
X = rectwin(L);
wrect = rectwin(L);
%apply window to signal
for n = 1:32
    yrect(n) = (80*dirac(n) - 80*sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16))*wrect(n);
end

N = L+1;
Xk = fftshift(fft(yrect,N));
figure(8);
plot(20*log10(abs(Xk)));

N = 10000;
Xk = fftshift(fft(yrect,N));
figure(9);
plot(20*log10(abs(Xk)));

L = 64;
%rectangular window
X = rectwin(L);
wrect = rectwin(L);
%apply window to signal
for n = 1:64
    yrect(n) = (80*dirac(n) - 80*sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16))*wrect(n);
end

N = L+1;
Xk = fftshift(fft(yrect,N));
figure(10);
plot(20*log10(abs(Xk)));

N = 10000;
Xk = fftshift(fft(yrect,N));
figure(11);
plot(20*log10(abs(Xk)));


L = 128;
%rectangular window
X = rectwin(L);
wrect = rectwin(L);
%apply window to signal
for n = 1:128
    yrect(n) = (80*dirac(n) - 80*sin(pi*n/2)/(pi*n) + cos(pi*n/32) + cos(pi*n/16))*wrect(n);
end

N = L+1;
Xk = fftshift(fft(yrect,N));
figure(12);
plot(20*log10(abs(Xk)));

N = 10000;
Xk = fftshift(fft(yrect,N));
figure(13);
plot(20*log10(abs(Xk)));

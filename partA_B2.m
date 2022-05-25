%Sapountzi Athanasia Despoina 2624
%Part A and Part B2

%% A1. FIR Filter
%methology for construction of FIR filter 
fsample = 16000;
As = 40;
Rp = 1;
fs = 0.15; 
fp = 0.05;
ws=0.3*pi; 
wp=0.1*pi;

%mininum filter length M from hamming's equation:
N = ceil( 4/(fs-fp) );

%discrete time fir-hamming-lowpass coefficients
coefficients = fir1(N, (wp+ws)/2, 'low', hamming(N+1));

%% A2. IIR Filter
%methology for construction of IIR filter Butterworth via bilinear transform
%First convert to analog
WP = 2*fsample*tan(wp/2);
WS = 2*fsample*tan(ws/2);

%Butterworth
[N,Wn] = buttord(WP, WS, Rp, As, 's');
[z, p] = butter(N, Wn, 'low', 's');

%bilinear transform
[num, den] = bilinear(z, p, fsample);


%% Impulse response
hold on;
figure (1);  
title('Impulse Response'); 
grid('ON');

ylabel('Amplitude'); 
xlabel('Samples');

%for FIR filter
[firY, firX] = impz(coefficients);
stem(firX, firY, 'm', 'filled');
hold on;
%for IIR filter
[iirY, IIR_x] = impz(num, den);
stem(IIR_x, iirY, 'b', 'filled');
hold on;
legend('FIR', 'IIR');
hold off;


%% Step response
figure (2);
hold on;
title('Step Response'); 
grid('ON');
xlabel('Samples');
ylabel('Amplitude');
%for FIR filter Step response
[FIR_y, FIR_x] = stepz(coefficients);
stem(FIR_x, FIR_y, 'm', 'filled');
hold on;
%for ΙIR filter Step response
[IIR_y, IIR_x] = stepz(num, den);
stem(IIR_x, IIR_y, 'b', 'filled');
hold on;
legend('FIR', 'IIR');
hold off;


%% Frequency response H(e^jw)
%frequency domain FIR
[FIR_h,FIR_w]=freqz(coefficients,1,256);
%amplitute-gain FIR
FIR_h_inDB = 20*log10(abs(FIR_h));
%frequency domain IIR
[IIR_h,IIR_w]=freqz(num,den,256);
%amplitute-gain IIR
IIR_h_inDB = 20*log10(abs(IIR_h));

figure(3); 
hold on; 
title('Frequency Response');
grid('ON'); 
ylabel('Gain');
xlabel('Normalised frequency');
%normalized frequency(Nyquist) FIR
plot(FIR_w/(pi), FIR_h_inDB, 'm');
hold on;
%Normalized frequency(Nyquist) IIR
plot(IIR_w/(pi), IIR_h_inDB,'b');
legend('FIR', 'IIR');
hold off;
%% Group delay
% -d/dw{ -tan^-1(  Im(H(e^jw))/Re(H(e^jw)) }
%by default there are 8192 samples
[FIR_gd, FIR_w] = grpdelay(coefficients);
[IIR_gd, IIR_w] = grpdelay(num, den);

figure(4); 
hold on; 
title('Group delay');
grid('ON');
xlabel('Normalised frequency');
plot(FIR_w/pi,FIR_gd, 'm');
hold on;
plot(IIR_w/pi, IIR_gd,'b');
legend('FIR', 'IIR');
hold off;

%% Zeroes/Poles
figure(5);
%only for IIR filter this time, because FIR does not have poles
[b, a] = zplane(z, p);
hold on; 
title('Zeros/Poles');
hold off;

%% B2.1
%for sample = 16KHz 
freqs = 0:10:5000;
% points per segment for 16KHz
window_1600 = hamming(1600); 
window_160 = hamming(160); 
% overlap with 80 points per segment (160/2=80)
%overlap is how many samples to include from the calculation of spectrum N-1 in spectrum N
overlap = 80; 

%read the record
[name_sample, fsample] = audioread('name_mono.wav');
%fsample=16000 afou h hxografhsh einai sta 16KHz.
%check the name_sample (if it is a vector)
whos mySpeech
%plot the recording
figure(10);
hold on; 
plot(name_sample); 
xlabel('samples=n');
hold off;

% hamming window L=1600
figure(11); 
hold on;
   
[S, F, T, P] = spectrogram(name_sample, window_1600, overlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time');
ylabel('Frequency(Hz)'); 
title('L=1600 Hamming window-5ms shift');  
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

% hamming window L=160
figure (12);
hold on;

[S, F, T, P] = spectrogram(name_sample, window_160, overlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time');
ylabel('Frequency(Hz)'); 
title('L=160 Hamming window-5ms shift');
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

%% B2.2 Sound spectrograms with lowpass(IIR)
%apply filter from part A
filtered_record = filter(num, den, name_sample);
audiowrite('filtered_name_16KHz.wav', filtered_record, fsample);


%Hamming window L=1600
figure (13); 
hold on; 
[S, F, T, P] = spectrogram(filtered_record, window_1600, overlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time'); 
ylabel('Frequency(Hz)'); 
title('L=1600 Hamming window-5ms shift');
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;

%Hamming window L=160
figure (14);
hold on; 
[S, F, T, P] = spectrogram(filtered_record, window_160, overlap, freqs, fsample, 'yaxis');
imagesc([0:0.1:20], F, 10*log10(abs(P)), [-160 -50] );
xlabel('Time'); 
ylabel('Frequency(Hz)'); 
title('L=160 Hamming window-5ms shift'); 
axis xy; axis tight; colormap(jet); view(0,90);
colorbar();
hold off;



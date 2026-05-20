clc;
clear;
close all;

%% PART I: ECG Signal Denoising - MIT-BIH Record

%% 1) Load ECG Record Manually

record = '100';

cd(['C:\Users\omars\Downloads\DSP_final_project_' record]); %rooh lel folder el fee el ecg

fs = 360;

fid = fopen([record '.dat'],'r'); %neftah el file 100.dat 3ashan ne2rah
if fid == -1
    error('Could not open ECG file.');
end

data = fread(fid, inf, 'uint8'); %ne2rah kol el data men el  file, inf read all, uint8 read as bytes
fclose(fid);

data = reshape(data,3,[])'; % han2asem el data le rows of 3 bytes we kol 3 bytes fee 2 ecg samples MIT-BIH bey5azen keda

s1 = data(:,1) + 256*bitand(data(:,2),15); % extract channel 1, bena5od byte 1, we el lower 4 bits beto3 byte 2 ne3mel 12-bit channel 1
s2 = data(:,3) + 256*bitshift(data(:,2),-4); % extract channel 2, bena5od byte 3, we el upper 4 bits beto3 byte 2 ne3mel 12-bit channel 2

s1(s1 >= 2048) = s1(s1 >= 2048) - 4096; % hawel men unsigned le signed
s2(s2 >= 2048) = s2(s2 >= 2048) - 4096; % same

gain = 200; % 200 counts= 1mv
adczero = 1024; % zero level

x = (s1 - adczero)/gain; % remove offset, convert to mv making real ECG

duration = 10; % 10s
N = duration * fs;
x = x(1:N);
t = (0:N-1)/fs;
figure;
plot(t,x);
grid on;
title(['Original ECG Record ' record]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');

%% 2) Frequency Spectrum Before Filtering

N = length(x);
X = fftshift(fft(x));
f = (-N/2:N/2-1)*(fs/N);

figure;
plot(f,abs(X)/N);
grid on;
title('Original ECG Frequency Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% 3) FIR WINDOW-BASED FILTER DESIGN

%remove baseline wander(low frequency<0.5hz)
fir_order_hp = 100; %kafy lel filter w mesh 3aly awy 3ashan may3melsh delay
b_fir_hp = fir1(fir_order_hp,0.5/(fs/2),'high'); %order, normalized cutoff(relative to fs/2), type
a_fir_hp = 1; %no feedback no den a=1

%remove powerline @50hz
fir_order_notch = 100; 
b_fir_notch = fir1(fir_order_notch,[49 51]/(fs/2),'stop');
a_fir_notch = 1;

%remove higher muscle noise 100-150hz as ecg from 0.5-100hz
fir_order_lp = 100;
b_fir_lp = fir1(fir_order_lp,100/(fs/2),'low');
a_fir_lp = 1;

%% 4) BUTTERWORTH IIR FILTER DESIGN

but_order = 4; %balanced

[b_but_hp,a_but_hp] = butter(but_order,0.5/(fs/2),'high');

[b_but_notch,a_but_notch] = butter(but_order,[49 51]/(fs/2),'stop');

[b_but_lp,a_but_lp] = butter(but_order,100/(fs/2),'low');

%% 5) CHEBYSHEV TYPE I IIR FILTER DESIGN

cheb_order = 4;
Rp = 0.5; %passband ripple small enough not to distort ECG shape significantly

[b_cheb_hp,a_cheb_hp] = cheby1(cheb_order,Rp,0.5/(fs/2),'high');

[b_cheb_notch,a_cheb_notch] = cheby1(cheb_order,Rp,[49 51]/(fs/2),'stop');

[b_cheb_lp,a_cheb_lp] = cheby1(cheb_order,Rp,100/(fs/2),'low');

%% 6) ANALYSIS OF ALL 9 FILTERS

%function to analyze filter plotting mag and phase response pole-zero plot
%impulse and step response
analyze_filter(b_fir_hp,a_fir_hp,'FIR High-pass Filter');
analyze_filter(b_fir_notch,a_fir_notch,'FIR Notch Filter');
analyze_filter(b_fir_lp,a_fir_lp,'FIR Low-pass Filter');

analyze_filter(b_but_hp,a_but_hp,'Butterworth High-pass Filter');
analyze_filter(b_but_notch,a_but_notch,'Butterworth Notch Filter');
analyze_filter(b_but_lp,a_but_lp,'Butterworth Low-pass Filter');

analyze_filter(b_cheb_hp,a_cheb_hp,'Chebyshev Type I High-pass Filter');
analyze_filter(b_cheb_notch,a_cheb_notch,'Chebyshev Type I Notch Filter');
analyze_filter(b_cheb_lp,a_cheb_lp,'Chebyshev Type I Low-pass Filter');

%% 7) APPLY FILTERS TO ECG

y_fir = filter(b_fir_hp,a_fir_hp,x);
y_fir = filter(b_fir_notch,a_fir_notch,y_fir);
y_fir = filter(b_fir_lp,a_fir_lp,y_fir); %cascaded filtering applies to remove all noise

y_but = filter(b_but_hp,a_but_hp,x);
y_but = filter(b_but_notch,a_but_notch,y_but);
y_but = filter(b_but_lp,a_but_lp,y_but);

y_cheb = filter(b_cheb_hp,a_cheb_hp,x);
y_cheb = filter(b_cheb_notch,a_cheb_notch,y_cheb);
y_cheb = filter(b_cheb_lp,a_cheb_lp,y_cheb);

%% 8) Time-Domain Comparison

figure;
subplot(4,1,1);
plot(t,x);
grid on;
title('Original ECG');
xlabel('Time (s)');
ylabel('mV');
ylim([-1 1.5]);

subplot(4,1,2);
plot(t,y_fir);
grid on;
title('Filtered ECG - FIR');
xlabel('Time (s)');
ylabel('mV');
ylim([-1 1.5]);

subplot(4,1,3);
plot(t,y_but);
grid on;
title('Filtered ECG - Butterworth');
xlabel('Time (s)');
ylabel('mV');
ylim([-1 1.5]);

subplot(4,1,4);
plot(t,y_cheb);
grid on;
title('Filtered ECG - Chebyshev Type I');
xlabel('Time (s)');
ylabel('mV');
ylim([-1 1.5]);

%% 9) Frequency Spectrum After Filtering

Y_fir = fftshift(fft(y_fir));
Y_but = fftshift(fft(y_but));
Y_cheb = fftshift(fft(y_cheb));

figure;

subplot(4,1,1);
plot(f,abs(X)/N);
grid on;
title('Original ECG Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-120 120]);
ylim([0 0.05]);

subplot(4,1,2);
plot(f,abs(Y_fir)/N);
grid on;
title('Filtered ECG Spectrum - FIR');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-120 120]);
ylim([0 0.05]);

subplot(4,1,3);
plot(f,abs(Y_but)/N);
grid on;
title('Filtered ECG Spectrum - Butterworth');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-120 120]);
ylim([0 0.05]);

subplot(4,1,4);
plot(f,abs(Y_cheb)/N);
grid on;
title('Filtered ECG Spectrum - Chebyshev Type I');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-120 120]);
ylim([0 0.05]);

%% 10) PSD Comparison

figure;

subplot(4,1,1);
pwelch(x,[],[],[],fs); %[] means Matlab default value so default window,window length # samples per segment,overlap,nfft and fs so axis in Hz instead of normalized
title('PSD of Original ECG');
ylabel('');

subplot(4,1,2);
pwelch(y_fir,[],[],[],fs);
title('PSD of Filtered ECG - FIR');
ylabel('');

subplot(4,1,3);
pwelch(y_but,[],[],[],fs);
title('PSD of Filtered ECG - Butterworth');
ylabel('');

subplot(4,1,4);
pwelch(y_cheb,[],[],[],fs);
title('PSD of Filtered ECG - Chebyshev Type I');
ylabel('');

%% 11) Spectrogram Comparison
%shows signal frequency content change over time
%color represents signal power yellow/bright=high energy, blue/dark=low
%energy
figure;

subplot(4,1,1);
spectrogram(x,128,120,128,fs,'yaxis'); %x → input signal, window length (samples per segment), overlap, nFFT points, fs  sampling frequency (gives Hz instead of normalized), 'yaxis'  frequency shown vertically, values best chosen for spectrogram
ylim([0 180]);
title('Spectrogram of Original ECG');

subplot(4,1,2);
spectrogram(y_fir,128,120,128,fs,'yaxis');
ylim([0 180]);
title('Spectrogram of Filtered ECG - FIR');

subplot(4,1,3);
spectrogram(y_but,128,120,128,fs,'yaxis');
ylim([0 180]);
title('Spectrogram of Filtered ECG - Butterworth');

subplot(4,1,4);
spectrogram(y_cheb,128,120,128,fs,'yaxis');
ylim([0 180]);
title('Spectrogram of Filtered ECG - Chebyshev Type I');

%% 12) SNR Improvement Estimate
%signal to noise ratio
noise_fir = x - y_fir; %original-filtered=noise
noise_but = x - y_but;
noise_cheb = x - y_cheb;

snr_fir = 10*log10(sum(x.^2)/sum(noise_fir.^2)); %squared to get power; ratio of signal power to noise power then converted to dB; higher snr better filtering
snr_but = 10*log10(sum(x.^2)/sum(noise_but.^2));
snr_cheb = 10*log10(sum(x.^2)/sum(noise_cheb.^2));

disp('Estimated SNR after FIR filtering in dB:');
disp(snr_fir);

disp('Estimated SNR after Butterworth filtering in dB:');
disp(snr_but);

disp('Estimated SNR after Chebyshev filtering in dB:');
disp(snr_cheb);

%% 13) Display All Filter Coefficients

disp(' FIR filter coefficients:');
disp('FIR High-pass b:');
disp(b_fir_hp);
disp('FIR Notch b:');
disp(b_fir_notch);
disp('FIR Low-pass b:');
disp(b_fir_lp);

disp('Butterworth filter coefficients:');
disp('Butterworth High-pass b:');
disp(b_but_hp);
disp('Butterworth High-pass a:');
disp(a_but_hp);
disp('Butterworth Notch b:');
disp(b_but_notch);
disp('Butterworth Notch a:');
disp(a_but_notch);
disp('Butterworth Low-pass b:');
disp(b_but_lp);
disp('Butterworth Low-pass a:');
disp(a_but_lp);

disp('Chebyshev filter coefficients:');
disp('Chebyshev High-pass b:');
disp(b_cheb_hp);
disp('Chebyshev High-pass a:');
disp(a_cheb_hp);
disp('Chebyshev Notch b:');
disp(b_cheb_notch);
disp('Chebyshev Notch a:');
disp(a_cheb_notch);
disp('Chebyshev Low-pass b:');
disp(b_cheb_lp);
disp('Chebyshev Low-pass a:');
disp(a_cheb_lp);

%% Filter Function

function analyze_filter(b,a,filter_name) %function defenition

imp = [1 zeros(1,200)]; %impulse signal for response
u = ones(1,200); %step signal for response

[H,w] = freqz(b,a,512); %frequency response 512 no of freq points(resolution)

h = filter(b,a,imp);
s = filter(b,a,u);

figure;

subplot(2,2,1);
plot(w/pi,20*log10(abs(H))); %mag response in dB standard and more visible, w/pi normalized
grid on;
title('Magnitude Response');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');

subplot(2,2,2);
plot(w/pi,angle(H));
grid on;
title('Phase Response');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Phase (rad)');

subplot(2,2,3);
stem(h);
grid on;
title('Impulse Response');
xlabel('n');
ylabel('Amplitude');

subplot(2,2,4);
plot(s);
grid on;
title('Step Response');
xlabel('n');
ylabel('Amplitude');

sgtitle(filter_name); %main title

figure;
zplane(b,a);
title([filter_name ' - Pole-Zero Plot']);

end
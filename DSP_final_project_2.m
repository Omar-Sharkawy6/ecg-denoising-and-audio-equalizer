clc;
clear;
close all;

[file, path] = uigetfile('*.wav', 'Select audio file');
[x, fs] = audioread(fullfile(path, file));
x = mean(x,2);

fprintf('Input sampling rate = %d Hz\n', fs);

mode = input('Enter mode (1 = Preset, 2 = Custom): ');

filter_structure = input('Filter structure (1 = FIR, 2 = IIR): ');

order = input('Enter filter order (recommended: FIR<=100, IIR<=10): ');

fs_out = input('Enter output sample rate: ');

if filter_structure == 1
    
    fprintf('\nFIR Window Types:\n');
    fprintf('1 = Hamming\n');
    fprintf('2 = Hanning\n');
    fprintf('3 = Blackman\n');
    
    fir_window = input('Select FIR window type: ');
    
else
    
    fprintf('\nIIR Filter Types:\n');
    fprintf('1 = Butterworth\n');
    fprintf('2 = Chebyshev Type I\n');
    fprintf('3 = Chebyshev Type II\n');
    
    iir_type = input('Select IIR filter type: ');
    
    if iir_type == 2
        Rp = input('Enter passband ripple Rp (dB): ');
    elseif iir_type == 3
        Rs = input('Enter stopband attenuation Rs (dB): ');
    end
    
end

if mode == 1 
    bands = [0 100;
             100 300;
             300 800;
             800 2000;
             2000 5000;
             5000 10000;
             10000 20000];
else
    
    nBands = input('Enter number of bands (5–10): ');
    bands = zeros(nBands,2);
    
    for i = 1:nBands
        bands(i,:) = input(['Enter band ', num2str(i), ' [f1 f2]: ']);
    end
    
end

nyq = fs / 2;

valid_idx = bands(:,1) < nyq;
bands = bands(valid_idx,:);

bands(:,2) = min(bands(:,2), nyq);

bands = bands(bands(:,2) > bands(:,1), :);

fprintf('Adjusted bands:\n');
disp(bands);

gains_db = zeros(size(bands,1),1);

for i = 1:length(gains_db)
    gains_db(i) = input(['Gain for band ', num2str(i), ' (dB): ']);
end

y_total = zeros(size(x));

for i = 1:size(bands,1)
    
    f1 = bands(i,1) / nyq;
    f2 = bands(i,2) / nyq;
    
    f1 = max(f1, 0.001);
    f2 = min(f2, 0.999);
    
    if f2 <= f1
        warning(['Skipping invalid band ', num2str(i)]);
        continue;
    end
    
    if filter_structure == 1
        
        switch fir_window
            case 1
                win = hamming(order+1);
                
            case 2
                win = hann(order+1);
                
            case 3
                win = blackman(order+1);
                
            otherwise
                error('Invalid FIR window selection');           
        end
        
        b = fir1(order, [f1 f2], 'bandpass', win);
        a = 1;
        
    else
        
        switch iir_type            
            case 1
                [b, a] = butter(order, [f1 f2], 'bandpass');
                
            case 2
                [b, a] = cheby1(order, Rp, [f1 f2], 'bandpass');
                
            case 3
                [b, a] = cheby2(order, Rs, [f1 f2], 'bandpass');
                
            otherwise
                error('Invalid IIR filter selection');                
        end        
    end
    
    y_band = filter(b, a, x);
    
    if any(isnan(y_band)) || any(isinf(y_band))
        warning(['Numerical issue in band ', num2str(i), ', skipping']);
        continue;
    end
    
    gain_linear = 10^(gains_db(i)/20);
    
    y_band = gain_linear * y_band;   
    y_total = y_total + y_band;
    
    figure;
    
    subplot(3,2,1);
    freqz(b,a,1024,fs);
    title(['Band ', num2str(i), ' Frequency Response']);    
    subplot(3,2,2);
    impz(b,a);
    title('Impulse');    
    subplot(3,2,3);
    stepz(b,a);
    title('Step');
    subplot(3,2,4);
    zplane(b,a);
    title('Pole-Zero');
    
end

if max(abs(y_total)) == 0
    error('Output signal is zero. Check bands or filter settings.');
end

y_total = y_total / max(abs(y_total));
y_out = resample(y_total, fs_out, fs);

t = (0:length(x)-1)/fs;

figure;
subplot(2,1,1);
plot(t, x);
title('Original Signal');
subplot(2,1,2);
plot(t, y_total);
title('Equalized Signal');

if any(isnan(y_total)) || any(isinf(y_total))
    error('y_total contains invalid values before PSD');
end

figure;

pwelch(x,[],[],[],fs);
hold on;

pwelch(y_total,[],[],[],fs);

legend('Original','Equalized');

title('Power Spectral Density');

figure;

subplot(2,1,1);
spectrogram(x,256,128,256,fs,'yaxis');
title('Original');
subplot(2,1,2);
spectrogram(y_total,256,128,256,fs,'yaxis');
title('Equalized');

disp('Playing original...');

sound(x, fs);

pause(length(x)/fs + 1);

disp('Playing equalized...');

sound(y_total, fs);

audiowrite('equalized_output.wav', y_out, fs_out);

disp('Done. Output saved as equalized_output.wav');
alpha = mod(205,3);
%{
%1st
N = 300
a = hann(N);
b= rectwin(N);
c = hamming(N);
figure;
plot(a)
title("hanning")
xlabel("sample index")
ylabel("magnitude")
figure;
plot(b)
title("rectangular")
xlabel("sample index")
ylabel("magnitude")
figure;
plot(c)
title("hamming")
xlabel("sample index")
ylabel("magnitude")
bl = blackman(N);
y = abs(fftshift(fft(bl,1024)))/N;
f = ((-length(y)/2:(length(y)/2-1)))/length(y);
b = 20*log(y);
figure;
plot(f,y)
title("spectrum N = 300")
xlabel("frequency")
ylabel("magnitude")
figure;
plot(f,b)
title("bode N = 300")
xlabel("frequency")
ylabel("magnitude")
%}
%{
%2nd
n =21;
wn = (pi/(alpha+1))*1/pi;
b = fir1(n-1,wn,"low",blackman(n));
c = fir1(n-1,wn,"low",rectwin(n));

figure;
impz(b)
title("blackman window FIR filter")
xlabel("sample index")
ylabel("magnitude")
figure;
freqz(b)
title("bode blackman FIR filter")
figure;
impz(c)
title("rectangular window FIR filter")
xlabel("sample index")
ylabel("magnitude")
hold on;
figure;
freqz(c)
title("bode rectangular FIR filter")
%}

%3rd
[x,Fs] = audioread("instru2.wav");
m = abs(fft(x));
F = (0:length(m)-1)*Fs/length(m);
figure;
plot(F,m)
title('instru2')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
fft1 = fft(x);
l1 = length(fft1);
mag1 = abs(fft1);
freq = (0:(l1-1))*(Fs/l1);
[~,peak] = max(mag1);
fund_f1 = freq(peak);
% Spectrogram parameters
window_length = 100; % Length of the Hamming window (samples)
overlap = 10;        % Overlap between consecutive windows (samples)
window_length_2 = 150;
fund_f = 775;
l = fund_f-40;
h = fund_f+40;
[A,B,C,D] = butter(8,[l h]/(Fs/2),"bandpass");
d1 = designfilt("bandpassiir",FilterOrder=8,HalfPowerFrequency1=l,HalfPowerFrequency2=h,SampleRate=Fs);
sos = ss2sos(A,B,C,D);
fvt = fvtool(sos);
legend(fvt,"butter")

filtered_audio = filtfilt(d1,x);
audiowrite('filtered_audio.wav', filtered_audio, Fs);
sound(filtered_audio, Fs);

% Create and plot the spectrogram

figure;
spectrogram(x, hamming(window_length), overlap);
title('Spectrogram of instru2');

figure;
[x1,Fs] = audioread("filtered_audio.wav");
spectrogram(x1, hamming(window_length), overlap);
title('Spectrogram of filtered output');

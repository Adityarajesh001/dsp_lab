clear
alpha = 1 + mod(205,3);

%1st 
t = -1:0.0001:1;
y1 = sin(2*pi.*alpha.*t);
y2 = cos((5*pi.*alpha.*t)+pi/4);
y3 = exp(-2*alpha.*t);
y4 = exp(-0.25*alpha.*t).*sin(20*pi*t);
figure(1);
xlabel("time");
ylabel("y1");
plot(t,y1);
figure(2);
xlabel("time");
ylabel("y1");
plot(t,y2);
figure(3);
xlabel("time");
ylabel("y1");
plot(t,y3);
figure(4);
xlabel("time");
ylabel("y1");
plot(t,y4);

%2nd
t = -5:0.001:5;
x = 46;
a = 1 + mod(x,3);
y = @(t)(exp(-a*t));
subplot(2,2,1);

plot(t,y(t));
title("1st");
subplot(2,2,2);

plot(t,y(-t));
title("2nd");
subplot(2,2,3);

plot(t,y(t-(1.5*a)));
title("3rd");
subplot(2,2,4);

plot(t,y(2*a*t));
title("4th");

%3nd
a = load('ECG_Data.txt');
plot(a);
b = load('RainFallIndia_Jan.txt');
histogram(b);
c = load('RainFallIndia_July.txt');
meam_jan = mean(b);
std_jan = std(b);
mean_july = mean(c);
std_july  = std(c);
fs = audioread('Track001.wav');
plot(fs)
xlabel('Time')
ylabel('Audio Signal')

%4th
F = 250.*a;
data = importdata('speech.wav');
[s,Fs] = audioread('speech.wav');
y1 = fn(s, F, Fs);
T = 1/Fs;       
L = 1500;
t = (0:L-1)*T;
y0 = fft(y1);
s0 = fft(s);
f = Fs*(0:(L/2))/L;
P2 = abs(y0/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Q2 = abs(s0/L);
Q1 = Q2(1:L/2+1);
Q1(2:end-1) = 2*Q1(2:end-1);
subplot(2,2,1);
plot(y1);
xlabel('time');
ylabel('y');
subplot(2,2,2);
plot(s);
xlabel('time');
ylabel('s');
subplot(2,2,3);
plot(f, P1,'blue', f, Q1,'red');%peaks with different magnitudes but occuring at same frequency
xlabel('frequency');
ylabel('Amplitude');
legend('y1','s')
function y = fn(s, F, Fs)
y = s.*cos(2.*pi.*F./Fs);
end


%5th
Fs = 44100;
alpha = 1 +mod(205,3);
t = 0.001:0.001:5;
y1 = sin(200*alpha*t);
y2 = sin(220*alpha*t);
y3 = y1+y2;
filename = 'exp1.wav';
audiowrite(filename,y1,Fs);
[y3,Fs] = audioread('exp2.wav');
sound(y3,Fs);
[y3,Fs]=audioread('exp2.wav');
t=linspace(0,length(y3)/fs,length(y3));
Nfft=16777216; %power of 2 and I put a huge number so there are many data point
f=linspace(0,fs,Nfft);
X1=abs(fft(y3,Nfft));
plot(f(1:Nfft/2),X1(1:Nfft/2))
xlabel('Frequency'); 
ylabel ('amp');
title ('FFT Spectrum');

% 6th
% Define frequencies for each note
freqs = [261.63, 293.66, 329.63, 349.23, 392.00, 440.00, 493.88, 523.25]; % Frequencies in Hz

% Define the duration of each note (in seconds)
note_duration = 0.5; % Duration of each note

% Create the time vector
fs = 44100; % Sampling frequency in Hz
t = 0:1/fs:note_duration-1/fs; % Time vector for one note

% Initialize the signal
signal = [];

% Generate the tones and append them to the signal
for i = 1:length(freqs)
    tone = sin(2 * pi * freqs(i) * t);
    signal = [signal, tone];
end

% Save the signal as a .wav file
filename = 'do_re_mi_fa_so_la_ti_do.wav';
audiowrite(filename, signal, fs);
y = audioread("do_re_mi_fa_so_la_ti_do.wav");
sound(y,fs)



%7th
%7th
load("ConvFile2.txt");
[y,Fs] = audioread('track002.wav');
C = conv2(ConvFile2,y);
filename = 'exp1_4.wav';
audiowrite(filename,C,Fs);
[C,Fs] = audioread('exp1_4.wav');

[y4,Fs]=audioread('exp1_4.wav');
t=linspace(0,length(y4)/Fs,length(y4));
Nfft=16777216; %power of 2 and I put a huge number so there are many data point
f=linspace(0,Fs,Nfft);
X1=abs(fft(y4,Nfft));
figure(1);
plot(f(1:Nfft/2),X1(1:Nfft/2))
xlabel('Frequency'); 
ylabel ('amp');
title ('FFT Spectrum');

[y,Fs]=audioread('Track002.wav');
t=linspace(0,length(y)/Fs,length(y));
Nfft=16777216; %power of 2 and I put a huge number so there are many data point
f=linspace(0,Fs,Nfft);
X1=abs(fft(y,Nfft));
figure(2);
plot(f(1:Nfft/2),X1(1:Nfft/2))
xlabel('Frequency'); 
ylabel ('amp');
title ('FFT Spectrum');
[y, Fs] = audioread('exp1_4.wav');
sound(y, Fs);

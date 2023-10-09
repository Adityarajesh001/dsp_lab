alpha = 2;
%{
Fs= 720;
Wp = 10/360;
Ws = 20/360;
[n,Wn] = buttord(Wp,Ws,2,40);
[b,a]=butter(n,Wn);
% tra
sys=tf(b,a,1/Fs)
figure;
zplane(b,a)
%text(real(z)+0.1,imag(z),"Zero")
%text(real(p)+0.1,imag(p),"Pole")
figure;
F = linspace(0,50,100)
[mag,phase] = bode(sys,2*pi*F)
plot(F,db(squeeze(mag)))
xlabel("freq")
ylabel("Mag(db)")
figure;
stepz(b,a);
title('Step Response');
figure;
impz(b,a);
title('Impulse Response');
%2nd
x = load("ECG_Data.txt")
y = filter(b,a, x)
t= (0: length(x)-1)/Fs;
figure;
plot(t,x)
hold on
plot(t,y)
legend('Input Data','Filtered Data')
ylabel("Amplitude")
xlabel("Time(s)")
%}

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
fund_f = freq(peak);
% Spectrogram parameters
window_length = 100; % Length of the Hamming window (samples)
overlap = 10;        % Overlap between consecutive windows (samples)
window_length_2 = 150;

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

%{
Fs= 720;
Wp_b = 10/360;
Ws_b = 20/360;
[n_b,Wn_b] = buttord(Wp_b,Ws_b,2,40);
[b_b,a_b]=butter(n_b,Wn_b);
sys_b=tf(b_b,a_b,1/Fs)
%figure;
%zplane(b_b,a_b)
%text(real(z)+0.1,imag(z),"Zero")
%text(real(p)+0.1,imag(p),"Pole")

Wp = 10/360;
Ws = 20/360;
[n,Wp] = cheb1ord(Wp,Ws,2,40)
[b,a]=cheby1(n,2,Wp);
% tra
sys=tf(b,a,1/Fs)

subplot(2,1,1)
zplane(b,a)
title("chebyshev")
subplot(2,1,2)
zplane(b_b,a_b)
title("butterworth")
%text(real(z)+0.1,imag(z),"Zero")
%text(real(p)+0.1,imag(p),"Pole")

figure;
F = linspace(0,50,100)
[mag,phase] = bode(sys,2*pi*F)
plot(F,db(squeeze(mag)))
legend("bode chebyshev")
xlabel("freq")
ylabel("Mag(db)")
hold on
[mag_b,phase_b] = bode(sys_b,2*pi*F);
plot(F,db(squeeze(mag_b)))
legend("bode chebyshev","bode butterworth")
hold off

figure;
impz(b,a);
hold on
impz(b_b,a_b);
legend("impulse chebyshev","impulse butterworth")
title('Impulse Response');
hold off
figure;
stepz(b,a);
hold on
legend("step chebyshev")
stepz(b_b,a_b);
legend("step chebyshev","step butterworth")
title('step Response');
hold off
%}
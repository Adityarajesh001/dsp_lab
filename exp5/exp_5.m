alpha = 1 + mod(205,3);
%{
t = 0:1/100:(10-1/100);
f = linspace(6,15,1000);
x = sin(2*pi.*f.*t);
figure(1)
plot(t,x)
title('chirp signal')
xlabel("time");
ylabel("AMPLITUDE");

m = abs(fft(x));
F = (0:length(m)-1)*100/length(m);
figure(2)
plot(F,m)
title('Magnitude')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");

% Spectrogram parameters
window_length = 100; % Length of the Hamming window (samples)
overlap = 10;        % Overlap between consecutive windows (samples)
window_length_2 = 150;

% Create and plot the spectrogram
figure(3);
spectrogram(x, hamming(window_length), overlap);
title('Spectrogram of Chirp Signal hamming 100');
figure(4);
spectrogram(x, hamming(window_length_2), overlap);
title('Spectrogram of Chirp Signal hamming 150');
figure(5);
spectrogram(x, hanning(window_length), overlap);
title('Spectrogram of Chirp Signal hanning');
figure(6);
spectrogram(x, blackman(window_length), overlap);
title('Spectrogram of Chirp Signal blackman');
%}
%{
[x,Fs] = audioread("instru2.wav");
m = abs(fft(x));
F = (0:length(m)-1)*Fs/length(m);
figure(1)
plot(F,m)
title('instru2')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");

% Spectrogram parameters
window_length = 100; % Length of the Hamming window (samples)
overlap = 10;        % Overlap between consecutive windows (samples)
window_length_2 = 150;

% Create and plot the spectrogram
figure(2);
spectrogram(x, hamming(window_length), overlap);
title('Spectrogram of instru2');


[x1,Fs1] = audioread("opera.wav");
m1 = abs(fft(x1));
F1 = (0:length(m1)-1)*Fs1/length(m1);
figure(3)
plot(F1,m1)
title('opera')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");

% Spectrogram parameters
window_length = 100; % Length of the Hamming window (samples)
overlap = 10;        % Overlap between consecutive windows (samples)

% Create and plot the spectrogram
figure(4);
spectrogram(x1, hamming(window_length), overlap);
title('Spectrogram of opera');
%}

[y,Fs] = audioread("name.wav");
figure(1)
plot(y)
title('name')
xlabel("time");
ylabel("AMPLITUDE");
m = abs(fft(y));
F = (0:length(m)-1)*Fs/length(m);
figure(2)
plot(F,m)
title('fft')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
window_length = 100; % Length of the Hamming window (samples)
overlap = 10;        % Overlap between consecutive windows (samples)
window_length_2 = 150;

% Create and plot the spectrogram
figure(3);
spectrogram(y, hamming(window_length), overlap);
title('Spectrogram of name');

%{
% Parameters
Fs = 10000;         % Sample rate (in Hz)
BitsPerSample = 16; % Bits per sample
Duration = 3;       % Duration of recording (in seconds)
FileName = 'name.wav'; % Output file name

% Create an audio recorder object
recObj = audiorecorder(Fs, BitsPerSample, 1);

% Record audio
disp('Start recording...');
record(recObj, Duration);
disp('Recording...');

% Wait for the recording to complete
pause(Duration);

% Stop recording
stop(recObj);
disp('Recording stopped.');

% Get the recorded audio data
audioData = getaudiodata(recObj);

% Save as a .wav file
audiowrite(FileName, audioData, Fs);

disp(['Audio saved as ' FileName]);
%}
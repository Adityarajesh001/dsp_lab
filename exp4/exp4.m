alpha = 1 + mod(205,3);


%1st

t = 0:1/120:2-1/120;
n = 15*alpha/120;
y = sin(30*2*pi*t);
m = abs(fft(y,120));
f = (0:length(m)-1)*120/length(m);
subplot(2,1,1)
plot(f,m)
title('Magnitude 120')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m1 = abs(fft(y,30));
f1 = (0:length(m1)-1)*120/length(m1);
subplot(2,1,2)
plot(f1,m1)
title('Magnitude 130')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
k = [1 2 3 4 5 6 7 8] ;
N1 = 120*k*n


%2nd
A = 140;
B = 146; 
t = 0:1/200:10-1/200;
xa = 0.1*sin(A*pi*t) + cos(B*pi*t);

m = abs(fft(xa,215));
f = (0:length(m)-1)*200/length(m);
figure(1)
plot(t,xa)
title('ORIGINAL SIGNAL')
xlabel("TIME");
ylabel("AMPLITUDE");
figure(2)

plot(f,m)
title('Magnitude 215')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m1 = abs(fft(xa,415));
f1 = (0:length(m1)-1)*200/length(m1);
figure(3)
plot(f1,m1)
title('Magnitude 415')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m2 = abs(fft(xa,1115));
f2 = (0:length(m2)-1)*200/length(m2);
figure(4)
plot(f2,m2)
title('Magnitude 1115')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m3 = abs(fft(xa,1515));
f3 = (0:length(m3)-1)*200/length(m3);
figure(5)
plot(f3,m3)
title('Magnitude 1515')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m4 = abs(fft(xa,1915));
f4 = (0:length(m4)-1)*200/length(m4);
figure(6)
plot(f4,m4)
title('Magnitude 1915')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");

%3rd
A = 140;
B = 146;
t = 0:1/200:10-1/200;
w = hann(2000);
x =(0.1*sin(A*pi*t) + cos(B*pi*t));
xa = (0.1*sin(A*pi*t) + cos(B*pi*t)).* w';
figure(1)
plot(t,xa)
title('EQUIVALENT SIGNAL')
xlabel("TIME");
ylabel("AMPLITUDE");
figure(2)
m = abs(fft(xa,215));
f = (0:length(m)-1)*200/length(m);
figure(2)
plot(f,m)
title('Magnitude 215')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m1 = abs(fft(xa,415));
f1 = (0:length(m1)-1)*200/length(m1);
figure(3)
plot(f1,m1)
title('Magnitude 415')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m2 = abs(fft(xa,1115));
f2 = (0:length(m2)-1)*200/length(m2);
figure(4)
plot(f2,m2)
title('Magnitude 1115')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m3 = abs(fft(xa,1515));
f3 = (0:length(m3)-1)*200/length(m3);
figure(5)
plot(f3,m3)
title('Magnitude 1515')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
m4 = abs(fft(xa,1915));
f4 = (0:length(m4)-1)*200/length(m4);
figure(6)
plot(f4,m4)
title('Magnitude 1915')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
figure(7)
plot(t,w)
title('Hanning')
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
figure(8)
plot(t,x)
title('ORIGINAL SIGNAL')
xlabel("TIME");
ylabel("AMPLITUDE");


%4th
load("Exp4Data2.txt")
x = Exp4Data2 ; 
m = abs(fft(x,10000));
f = (0:length(m)-1)*1/length(m);
figure(1)
plot(f,m)
title("WITHOUT HAMMING")
xlabel("FREQUENCY");
ylabel("AMPLITUDE");

w = hamming(500);
x1 = Exp4Data2.* w' ; 
m1 = abs(fft(x1,10000));
f1 = (0:length(m1)-1)*1/length(m1);
figure(2)
plot(f1,m1)
title("WITH HAMMING")
xlabel("FREQUENCY");
ylabel("AMPLITUDE");
figure(3)
plot(w)
title("HAMMING")
xlabel("POINTS");
ylabel("AMPLITUDE");


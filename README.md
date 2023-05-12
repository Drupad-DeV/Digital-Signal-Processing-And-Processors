# signal
## labsheet 5 
## Overlap Add:
Determine the output of an LTI system with impulse response h[n]={1,-1,1,2} for the input 
x[n]={1,2,-1,3,2,1-2,1,0,-1,-2,3,-2,1,1,2,2,1} using overlap add method
<pre>
  clear all
clc
x=[1 2 -1 3 2 1 -2 1 0 -1 -2 3 -2 1 1 2 2 1]
h=[1 -1 1 2]
xn=length(x);
hn=length(h);

if rem(xn,hn) ~= 0
    x=[x zeros(1,hn-rem(xn,hn))];
end

xn=length(x);
i=xn/hn;
a=zeros(1,xn+hn-1);
b=a;
H=fft(h,hn+hn-1);
for j=1:i
    X=fft(x((j-1)*hn+1:(j)*hn),hn+hn-1);
    z=H.*X;
    z=ifft(z,hn+hn-1);
    b=zeros(1,xn+hn-1);
    nz=length(z);
    b((j-1)*hn+1:nz+(j-1)*hn)=z;
    a=a+b;
end   
round(a);
disp (a);
</pre>

## Overlap Save:
Determine the output of an LTI system with impulse response h[n]={1,-1,1,2} for the input 
x[n]={1,2,-1,3,2,1-2,1,0,-1,-2,3,-2,1,1,2,2,1} using overlap save method.
<pre>
clear all
clc
x=[1 2 -1 3 2 1 -2 1 0 -1 -2 3 -2 1 1 2 2 1]
h=[1 -1 1 2]

xn=length(x);
hn=length(h);
N = xn+hn-1;
h1 = [h zeros(1,N-xn)];
H = fft(h1);
n = length(h1);
a = zeros(1,N);
x1 = [zeros(1,n-hn) x zeros(1,n)];
for i = 1:hn:N
    y = x1(i:i+(2*(n-hn)));
    y1 = fft(y);
    y2 = y1.*H; 
    y3 = round(ifft(y2));
    a(i:(i+n-hn)) =y3(hn:n);
end
disp(a)
</pre>

## Labsheet 6 Spectrum Using FFT
  2.Generate a sine wave with fm=200Hz and Fs=1000Hz. Plot the spectrum of this sine wave using FFT. (x axis 
should be frequency in Hz) 
<pre>
t = 0:1/1000:1-1/1000;  
fm = 200;
figure;
x = sin(2*pi*fm*t);
plot(t, x);
title("input Sine wave");
xlabel('Time');
ylabel('Magnitude');

% spectrum of the sine wave
N = length(x);
X = fft(x);
f = (0:N-1)*(1000/N);  
figure;
plot(f,abs(X));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum of Sine Wave');

for i =1:length(X)/2
 if (X(i)==max(X(1:(1000/2)+1)))
 dfrequency=i-1;
 break
 end
end
disp(dfrequency);
</pre>

1. 8 Point DFT into 4 point and 2 point
<pre>
clc
clear
x = [1 2 3 4 5 6 7 8];
x1(1) = x(1) + x(5);
x1(2) = x(1) - x(5);
x1(3) = x(3) + x(7);
x1(4) = x(3) - x(7);
x1(5) = x(2) + x(6);
x1(6) = x(2) - x(6);
x1(7) = x(4) + x(8);
x1(8) = x(4) - x(8);
x2(1) = x1(3) + x1(1);
x2(2) = x1(4) + x1(2).*(-j);
x2(3) = x1(1) - x1(3);
x2(4) = x1(2) - x1(4).*(-j);
x2(5) = x1(7) + x1(5);
x2(6) = x1(8) + x1(6).*(-j);
x2(7) = x1(5) - x1(7);
x2(8) = x1(6) - x1(8).*(-j);
%final output
y(1) = x2(5) + x2(1);
y(2) = x2(6)*(0.707 -0.707j) + x2(2);
y(3) = x2(7)*-j + x2(3);
y(4) = x2(8)*(-0.707 -0.707j) + x2(4);
y(5) = x2(1) - x2(5);
y(6) = x2(2) - x2(6)*(0.707 -0.707j);
y(7) = x2(3) - x2(7)*-j;
y(8) = x2(4) - x2(8)*(-0.707-0.707j);
disp(y);
</pre>

3. Audio SIgnal
<pre>
clear all
close all
clc
fm = 200; 
n = 0:1/1000:1; 
[x ,fs]=audioread('sample.au') ;
figure;
plot(1:length(x), x);
title('Input Signal')
xlabel('Time (s)');
ylabel('Amplitude');
N = length(x); 
X = abs(fft(x))
figure;
plot(0:N-1,abs(X));
title('Spectrum');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
</pre>

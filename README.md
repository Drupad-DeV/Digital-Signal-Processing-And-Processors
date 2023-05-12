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





##  Labshet 3: DFT Coefficients
Consider a Sequence 𝑥[𝑛] = ∑ 𝛿(𝑛 − 2𝑘)
a. Sketch x[n] for 0 ≤ 𝑛 ≤ 100
b. Find the DFS coefficients. 
c. Plot its Magnitude and Phase spectrum using Parseval’s relation. 
2. Consider 𝑥[𝑛] = cos(pi/3 * n) + sin (pi/4 * n)
. Determine its DFS coefficients. 
(2 more questions will be given on the spot for evaluation) 
<pre>
%Signal Generation
n=0:2:100;
x=ones(1,51);
stem(n,x)
xlabel("Time")
ylabel("Amplitude")
title("Signal Generation")
%DFS coefficient
sum1=0;
N=2;
for k=0:N-1
  for n=0:N-1
       sum1=sum1+(1/N*x(n+1)*exp(-1i*k*pi*n));
  end
    D(k+1)=sum1;
    sum1=0;
end
disp(D)
%Parseval's Relation
c=abs(D(k));
figure
stem(k,c)
xlabel("K")
ylabel("Dk")
title("Parseval's Relation")
</pre>

2. Consider 𝑥[𝑛] = cos (pi/3 * n) + sin (pi/4 * n)
 Determine its DFS coefficients. 
<pre>
%DFS Coefficients for Cos(pi/3*n)+Sin(pi/4*n)
x=0;
N=24;
for k=0:N-1
  for n=0:N-1
       x=x+(1/N*(sin(n*pi/4)+cos(n*pi/3)*exp(-1i*k*pi*n/12)));
  end
    D(k+1)=x;
    x=0;
end
disp(D)
stem(D)
</pre>

## Labsheet 4 DFT
1. Sketch the signal 𝒙[𝒏] = 𝒔𝒊𝒏ቀ𝟐 ∗ 𝒑𝒊 ∗
[𝟎:𝟏𝟗] / 
𝟐𝟎 
2. Find the DFT of x[n]. 
3. Shift x[n] by circularly by 2 units and name it as 𝒙𝟏
[𝒏]. Plot 𝒙𝟏
[𝒏]. 
4. Find DFT of 𝒙𝟏
[𝒏]. Correlate the relationship between X(k) and 𝑿𝟏(𝒌)

<pre>
N = 20;
n = 0:19;
w = 2*pi/N;
D = zeros(1,20);
D1 = [];
x = sin(2*pi*n/N);
subplot(3,2,1)
stem(n,x)
xlabel('time')
ylabel('Amplitude')
title('sin(2*pi*n/N)')
for k = n
 D(k+1) = sum(x.*exp(-1i*k*2*pi*n/N));
end
subplot(3,2,2)
stem(n,abs(D))
xlabel('samples')
ylabel('Amplitude')
title('DFT of x[n]')
x1 = x(mod(n-2,N)+1);
subplot(3,2,3)
stem(n,x1)
xlabel('time')
ylabel('Amplitude')
title('x1[n]')
for k = n
 D1(k+1) = sum(x1.*exp(-1i*k*2*pi*n/N));
end
subplot(3,2,4)
stem(n,abs(D1))
xlabel('samples')
ylabel('Amplitude')
title('DFT of x1[n]')
subplot(3,2,5)
stem(n,angle(D))
xlabel("time")
ylabel("Amplitude")
title("Phase Angle of x[k]")
subplot(3,2,6)
stem(n,angle(D1))
xlabel("time")
ylabel("Amplitude")
title("Phase Angle of x1[k]")
</pre>
2. The DFT of the 5-point signal x[n] is given by 𝑿[𝒌] = {𝟓, 𝟔, 𝟏, 𝟐, 𝟗}. Another signal is defined by 
𝒙𝟏
[𝒏] = 𝑾𝟓
^𝟐𝒏𝒙[𝒏]; 𝟎 ≤ 𝒏 ≤ 𝟒. Determine X1[k] using suitable property. Plot x[n], Real part of 
X[k], Imaginary part of X[k] in one figure. Plot x1[n], Real part of X1[k], Imaginary part of X1[k] in 
another figure. 

<pre>
n = 0:4;
X = [5, 6, 1, 2, 9];
N = 5;
X1 = X(mod(n-2,N)+1);
x = [];
x1 = [];
for k = n
   x(k+1) = sum(X.*exp(1i*k*2*pi*n/N))*1/N;
end
for k = n
   x1(k+1) = sum(X1.*exp(1i*k*2*pi*n/N))*1/N;
end
figure(1)
subplot(3,1,1)
stem(n, abs(X), 'r');
xlabel("n")
ylabel("values")
title("Magnitude of X")
subplot(3,1,2)
stem(n, angle(X), 'g');
xlabel("n")
ylabel("values")
title("Phase of X")
subplot(3,1,3)
stem(n, abs(x), 'b');
xlabel("n")
 
ylabel("values")
title("Magnitude of x")
figure(2)
subplot(3,1,1)
stem(n, abs(X1), 'r');
xlabel("n")
ylabel("values")
title("Magnitude of X1")
subplot(3,1,2)
stem(n, angle(X1), 'g');
xlabel("n")
ylabel("values")
title("Phase of X1")
subplot(3,1,3)
stem(n, abs(x1), 'b');
xlabel("n")
ylabel("values")
title("Magnitude of x1")
figure (3)
y1 = x.*exp(-1i*4*pi*n/N);
figure(3);
subplot(2,1,1);
stem(n, abs(y1));
title('Magnitude of y1');
xlabel('n');
ylabel('|y1|');
subplot(2,1,2);
stem(n, angle(y1));
title('Phase of y1');
xlabel('n');
ylabel('∠y1');
</pre>

3. Let x[n] be the finite sequence and its DFT is 𝑿[𝒌] = {𝟎, 𝟏 + 𝒋, 𝟏, 𝟏 − 𝒋}.Using the properties, find 
DFT of the following sequences. 
a. 𝒙𝟏
[𝒏] = 𝒆
𝒋𝝅
𝟐
𝒏𝒙[𝒏]
b. 𝒙𝟐
[𝒏] = 𝒄𝒐𝒔(𝒏𝝅
𝟐
)𝒙[𝒏]
c. 𝒙𝟑
[𝒏] = 𝒙(𝒏 − 𝟏)
<pre>
n = 0:3;
X = [0,1+1i,1,1-1i];
N = 4;
for k = n
    x(k+1) = sum(X.*exp(1i*k*2*pi*n/N))*1/N;
end


X1 = X(mod(n-1,N)+1);
for k = n
 x1(k+1) = sum(X1.*exp(1i*k*2*pi*n/N))*1/N;
end


X2 = X(mod(n+1,N)+1);
for k = n
 x2(k+1) = sum(X2.*exp(1i*k*2*pi*n/N))*1/N;
end


X3 = X.*exp(-1i*n*pi*1/N);
for k = n
 x3(k+1) = sum(X3.*exp(1i*k*2*pi*n/N))*1/N;
end


x1_n = x.*exp(1i*pi*n/2);
x2_n = x.*cos(n*pi/2);
x3_n = x(mod(n-1,4)+1);
abs(X)
abs(X1)
abs(X2)
abs(X3)
</pre>  

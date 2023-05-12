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


function b=BayesPSDb2(d,fs,wmean)
N=length(d);
t=0:1/fs:(N-1)/fs;

A=(d.*t)*(sin(wmean*t))';
B=(d.*t)*(cos(wmean*t))';
C=(d*(cos(wmean*t))')*((d.*t.*t)*(cos(wmean*t))');
D=(d*(sin(wmean*t))')*((d.*t.*t)*(sin(wmean*t))');
b=-(A^2+B^2-C-D)*2/N;
close all; clear; clc;

cd '/Users/zymantas/TVPVARPkg/test/TestOls'

T=100;
N=7;
x=randn(T,N);
beta=1:N;
y=x*beta' + randn(T,1);

data=[y x];
dlmwrite('data2test.txt',data,'precision',35)

x\y

[B,BINT,R,RINT,STATS] = regress(y,x);

s = regstats(y,x,'linear');

STATS=[s.rsquare s.mse];

dlmwrite('B2check.txt',s.beta,'precision',36);
dlmwrite('stats2check.txt',STATS,'precision',36);
dlmwrite('covb2check.txt',s.covb,'precision',36);




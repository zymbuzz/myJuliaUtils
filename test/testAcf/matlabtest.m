close all;clear; clc;

cd '/Users/zymantas/TVPVARPkg/test/testAcf/'

k=10;

y=[]
a=[]
b=[];
for i=1:k
x = randn(1000,1);         % 1000 Gaussian deviates ~ N(0,1)
z = filter([1 -1 1],1,x)
y = [y z];  % Create an MA(2) process
a= [a autocorr(z,'NumLags',40)];      % Inspect the ACF with 95% confidence
b=[b myautocorr(z,40)];
end

a==b
all(round(a-b,14)==0)

dlmwrite('y.txt',y,'delimiter',',','precision',25)
dlmwrite('acf.txt',b,'delimiter',',','precision',25)

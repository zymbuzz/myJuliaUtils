close all; clear; clc;
cd '/Users/zymantas/TVPVARPkg/test/TestWish/'

K=20;
T=80;

if T>81
    error('note that function for wishrnd is different if T>81') 
end

rng(50)
aa=randn(K,K);
aa=aa'*aa;
dlmwrite('WW.txt', aa ,'delimiter', ',' , 'precision', 15)
dlmwrite('T.txt', T ,'delimiter', ',' ,'precision', 15)

rng(50)
draw=randn(T,K);
dlmwrite("draw.txt",draw,'delimiter',',','precision',15)

rng(50)
matlabwish=wishrnd(aa,T);
dlmwrite("matlabwish.txt",matlabwish,'delimiter',',','precision',15)

rng(50)
matlabinvwish=iwishrnd(aa,T);
dlmwrite("matlabinvwish.txt",matlabinvwish,'delimiter',',','precision',15)


rng(50)
sigma=aa
df=T
[n,m] = size(sigma);
[d,p] = cholcov(sigma,0)
x = randn(df,n);

[~,R] = qr(x,0);
T = d' / R;

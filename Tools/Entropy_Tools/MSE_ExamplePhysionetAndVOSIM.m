%Example
N=10000;
% Simulate Stochastic 
noise=randn(N,1);
%Simulate determistic system with noise-like 2nd order statistics
nlinear=zeros(N,1);nlinear(1)=0.2;u=4;
for n=2:N;nlinear(n)=u*nlinear(n-1)*(1-nlinear(n-1));end

% MSE paramneters
maxScale=10;
r = 0.15; % default in Physione
m = 2;    % default in Physione

% VOSIM implementation 
[entropyNoiseVOSIM]=ComputeMultiscaleEntropy(noise,m,r,maxScale);
[entropyDetermVOSIM]=ComputeMultiscaleEntropy(nlinear,m,r,maxScale);
figure(1)
subplot(3,1,1);
plot(noise(1:1000));hold on;grid on;plot(nlinear(1:1000),'r');legend('Stochastic','Deterministic')
subplot(3,1,2);
plot(entropyNoiseVOSIM);hold on;grid on;plot(entropyDetermVOSIM,'r');legend('Stochastic','Deterministic')


% Using Physionet
[entropyNoise,scale1]=msentropy(noise,[],[],[],[],[],[],[],maxScale);
[entropyDeterm,scale2]=msentropy(nlinear,[],[],[],[],[],[],[],maxScale);
subplot(3,1,3);
plot(scale1,entropyNoise);hold on;grid on;plot(scale2,entropyDeterm,'r');legend('Stochastic','Deterministic')



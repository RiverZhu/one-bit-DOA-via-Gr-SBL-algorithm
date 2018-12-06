% Onebit DOA via the generalized SBL algorithm
% This code is written by Jiang Zhu and Xiangming Meng. If you have any
% problems, please feel free to contact Jiang Zhu via jiangzhu16@zju.edu.cn
% Paper: X. Meng and J. Zhu, A generalized sparse Bayesian learning
% algorithm for one-bit DOA estimation,  IEEE Communications Letters, 
% vol. 22, no. 7, pp. 1414-1417, 2018.

clc;  clear;   close all;
rng(6)      % 2,3,4,5,6 work
N = 361;                % grid size
maxit_outer = 400;
supp=[181-6  181+4 181+75*2];
K = length(supp);
x_dB = [12;22;20];  % amplitudes
% Bearing grid
theta = (-90:180/(N-1):90);
theta_r = theta*pi/180;
u = sin(theta_r);
d = 1/2;                % intersensor spacing
SNRdB = 40;             % [40]
M = 256;
L = 1;
c_sign = @(cpl_num)sign(real(cpl_num))+1j*sign(imag(cpl_num));
wvar0 = 1;

q = 0:1:(M-1);          % sensor numbering
xq = (q-(M-1)/2)*d;     % sensor locations
A = exp(-1i*2*pi*xq'*u)/sqrt(N); % M*N
x_amp = 10.^(x_dB/20);
x_amp = x_amp*ones(1,L);
X = zeros(N,L);
X(supp,:) = x_amp.*exp(1j*2*pi*rand(K,L));
if(L>1)
     wvar = ((norm(A*X,'fro'))^2/M/L)*10^(-SNRdB/10); 
else
     wvar = (norm(A*X))^2/M*10^(-SNRdB/10); 
end
% noise generation
w = sqrt(wvar/2)*randn(M,L)+1i*sqrt(wvar/2)*randn(M,L);
Y = c_sign(A*X+w);
% noise variance match
[theta_uninfor, NMSE_SBL,X_debiased] = onebitdoa_uninfor_iter( N, L, M, K, X, Y, wvar, maxit_outer );
% noise variance mismatch
[theta_uninfor_mismatch, NMSE_SBL_mismatch,X_debiased1] = onebitdoa_uninfor_iter( N, L, M, K, X, Y, wvar, maxit_outer );

theta_uninfor
theta_uninfor_mismatch
figure(1)
semilogx(1:maxit_outer,NMSE_SBL,'-b+',...
    1:maxit_outer,NMSE_SBL_mismatch,'-r<')
xlim([1 maxit_outer])
legend('matched','mismatched')
xlabel('iterations')
ylabel('debiased NMSE (dB)')

figure(2)
stem(abs(X_debiased),'ro')
hold on
stem(abs(X_debiased1),'b*')
legend('matched','mismatched')



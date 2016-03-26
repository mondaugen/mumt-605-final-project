function [ap,mp,thp,op,psp] = gdmr(x,Fs)
% Refined generalized derivative method as described by Hamilton (see
% references).
% x is the signal to analyse
% Fs is the sampling rate.
% By NICHOLAS ESTERER
M = 1023;
N = length(x);
n = [0:N-1];
t = n / Fs;
dx = derivative_filter(x,M,Fs);
ddx = derivative_filter(dx,M,Fs);
w = horzcat(hanning(N-1)',[0]); % assumes length of x is even
tw = w .* t;
o = linspace(0,2*pi,N+1)(1:end-1);
X = fft(x .* w);
DX = fft(dx .* w);
% find maximum magnitude value and index
[ap0,idx] = max(abs(X));
% get bin frequency
%op0 = (idx - 1) / N * 2 * pi;
% calculate frequency estimate for bin with greatest amplitude
op = imag(DX(idx) ./ X(idx));
% get stft of dx at estimated frequency
DXo = stcft(dx,op,w,1/Fs);
% get stft of x at estimated frequency
Xo = stcft(x,op,w,1/Fs);
% get stft of ddx at estimated frequency
DDXo = stcft(ddx,op,w,1/Fs);
% get stft of twx
TXo = stcft(x,op,tw,1/Fs);
% get stft of twdx
DTXo = stcft(dx,op,tw,1/Fs);
% get refinement factor
up = real(DTXo/Xo) - real(TXo * DXo / (Xo * Xo));
% calculate amplitude modulation estimate
mp = real(DXo / Xo);
% calculate frequency slope estimate
psp = (imag(DDXo / Xo) - 2*mp*op) / (1 + up);
Go = gamma_init_est(0,mp,psp,w,1/Fs,N);
ap = abs(Xo/Go);
thp = arg(Xo/Go);
endfunction

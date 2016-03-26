function [Xo] = stcft(x,o,w,Ts)
% STCFT short time continuous fourier transform
% evaluates sum_{n=-N/2}^{N/2} x[n+N/2] * w[n] * exp(-j*o*n)
% x is signal, o is frequency to investigate, w is window and Ts is sample
% period
% by NICHOLAS ESTERER
N = length(x);
if (mod(N,2) == 1)
    n = [-1*(N-1)/2:(N-1)/2]; % N odd
else
    n = [-1*N/2:N/2-1]; % N even
end
Xo = 2 / N * sum(x .* w .* exp(-j * o * Ts .* n));
endfunction

function [Xo] = gamma_init_est(o,u0,ps0,w,Ts,N)
% Evaluation of the gamma function as described by Depalle and Marchand
% evaluates sum_{n=-N/2}^{N/2} w[n] * exp(u0*t + j(o*t + ps0/2*(t*t)))
% x is signal, o is frequency to investigate, w is window and Ts is sample
% period
% by NICHOLAS ESTERER
if (mod(N,2) == 1)
    n = [-1*(N-1)/2:(N-1)/2]; % N odd
else
    n = [-1*N/2:N/2-1]; % N even
end
t = n * Ts;
Xo = 2 / N * sum(w .* ...
    exp(u0*t + j * ( o * t + ps0 / 2 * (t .^ 2))));
endfunction

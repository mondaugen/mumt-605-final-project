function [y] = derivative_filter(x,N,Fs)
% DERIVATIVE_FILTER finds an approximation of the derivative of x by convolving
% with a N-point differentiating filter (see Depalle and Marchand). Signal x was
% sampled at a rate Fs in Hz.
% By NICHOLAS ESTERER
n = [1:(N-1)/2]; % assume N odd
n_rev = fliplr(n) * -1;
h = [0];
h = horzcat(Fs * ((-1) .^ n_rev) ./ n_rev,h);
h = horzcat(h,Fs * ((-1) .^ n) ./ n);
h = h .* hanning(N)';
y = conv(h,x)((N-1)/2 + 1 : end - (N-1)/2); % remove head and tail of convolution
endfunction

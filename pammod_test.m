clear
clc
%% paramter definitions
N = 1e5;   % number of input bits
M = 32;      % modulatioin order
n = log2(M);    % number of bits per symbol
n_sym = floor(N/n);
N_bit = n_sym * n;
data = randi([0, 1], 1, N_bit);

%% transmitter
data_reshape = reshape(data, [n_sym, n]);

data_dec = bi2de(data_reshape);
signal = pammod(data_dec, M);

%% receiver
SNR = 0:15;

for snr_index = 1:length(SNR)
    signal_noisy = awgn(signal, SNR(snr_index), 'measured');
    sym_demod = pamdemod(signal_noisy, M);
    data_bin = de2bi(sym_demod);
    data_out = reshape(data_bin, [1, N_bit]);
    [~, ber(snr_index)] = biterr(data, data_out);
end

semilogy(SNR, ber, 'linewidth', 2);
grid on
xlabel('SNR (dB)')
ylabel('BER')
title('BER vs. SNR (PAM)')
hold on

clear all;
close all;
clc;

% LDPC Investigation

%% LDPC Code Simulation: Constants
SNRs = 4:.5:6;  % SNR range
R = 5/6;             % Code rate
iter = 10^2;         % Number of codewords per SNR                                                                                     
max_iter = 10;       % Maximum number of iteration for the BP algorithm

% Initialize Array
BER_BP = zeros([1 length(SNRs)]);
BER_BP_custom = zeros([1 length(SNRs)]);
nErrors_BP = zeros([1 length(SNRs)]);
nErrors_BP_custom = zeros([1 length(SNRs)]);

P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
blockSize = 27;
H = ldpcQuasiCyclicMatrix(blockSize,P);
% blockSize = 3;
% p = [0 -1 1 2; 2 1 -1 0];
% H = ldpcQuasiCyclicMatrix(blockSize,p);

% Create LDPC encoder and decoder configuration objects
cfgLDPCEnc = ldpcEncoderConfig(H);  
cfgLDPCDec = ldpcDecoderConfig(H);

%% Belief Propagation Algorithm
for idx_SNR = 1:length(SNRs)
    for idx_iter = 1:iter
        SNR = SNRs(idx_SNR);
        SNRlin = 10^(SNR/10); % Signal-to-Noise ratio for the current iteration
    
        % Generate bitstream
        info_bits = randi([0 1],cfgLDPCEnc.NumInformationBits,1);
    
        % Encode bitstream
        bits_enc = ldpcEncode(info_bits,cfgLDPCEnc);
    
        % BPSK modulation: s=-1 -> bit=0 and s=1 -> bit=1
        X = 2*bits_enc-1;
    
        % Transmission over bi-AWGN channel
        noise = randn(size(X));     % Noise generation    
        % Scale the noise power
        noise_scaled = sqrt(1/(SNRlin)) * noise;
        
        % Add AWGN
        Y = X + noise_scaled;

        % BPSK Demodulation
        X_hat = Y >= 0;

        if(sum(X_hat ~= bits_enc) > 0)
            % LLR Calculation
            sigma2 = 1/SNRlin;
            L = -2*Y/sigma2;
    
            % Matlab's BP algorithm
            bits_hat = ldpcDecode(L, cfgLDPCDec, max_iter); 

            % Custom LDPC BP algorithm
            bits_hat_c = ldpc_bp_function(L, H, max_iter);
        end

        nErrors_BP(idx_SNR) = nErrors_BP(idx_SNR) + biterr(bits_hat, info_bits);
        nErrors_BP_custom(idx_SNR) = nErrors_BP_custom(idx_SNR) + biterr(bits_hat_c, bits_enc);
    end
end

BER_BP = nErrors_BP./(iter * cfgLDPCEnc.NumInformationBits);
BER_BP_custom = nErrors_BP_custom./(iter * 648);

%% Plotting
figure;
semilogy(SNRs, BER_BP, '.-');
hold on;
semilogy(SNRs, BER_BP_custom, '.-');
set(gca,'XTick', SNRs);
xlabel('SNR [dB]')
ylabel('BER')
ylim([10^-4 10^0])
grid on;
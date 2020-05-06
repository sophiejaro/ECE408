% Sophie Jaro
% ECE408: Wireless Communications
% Part 1: 2x2 MIMO Link

clear all;
clc;

% Part 1a: 2x2 MIMO System with Zero-Forcing

% This program simulates a 2x2 MIMO system with 1000 symbols. The signal is
% modulated with BPSK, 4 QAM, 16 QAM, and 64 QAM to compare BER and data
% rates. The signal is passed through three different flat fading Rayleigh 
% channels, generated with max doppler shifts of 10, 100 and 1000 Hz. The
% signal is equalized with zero forcing to minimize ISI. Assuming channel
% state information at the reciever (CSIR), the inverse of the channel gain 
% % matrix is applied to the recieved signal. 
% The simulation is repeated 10,000 times, and the trials are averaged. 
% My goal was to obtain a bit error rate on the order of 1e-4 
% for an SNR of 20 and doppler frequency shift of 10 Hz. A pattern of BER
% increasing with the order of QAM and the max doppler frequency shift of
% the channel was observed, which makes sense since higher orders of QAM increase 
% data rate (the bits transmitted per symbol) but also increase the bit
% error rate. Also, the larger max doppler frequency shift causing a worse
% BER because the channel introduces more distortion to the signal. 

% Parameters
N = 2^0; % number of frequency domain points for Rayleigh Channel generation
nSym = 1000; % number of symbols per packet
numIter = 10000; % number of iterations

EbNo = 0:25;
ebno_index = 21;
SNR = 20;

fm = [10 100 1000]; % max doppler frequency shift (in Hz)
M = [2 4 16 64];

for i = 1:length(fm)% LOOP 1 for Max Doppler Frequency Shift
    for j = 1:length(M)% LOOP 2 for M-ary number
        for k = 1:numIter% LOOP 3 for Number of Iterations

            % Get Rayleigh Channel Coefficient Matrix
            h11 = raychan(fm(i),N);
            h12 = raychan(fm(i),N);
            h21 = raychan(fm(i),N);
            h22 = raychan(fm(i),N);
            H = [h11 h12; h21 h22];
            
            Mary = M(j);
            
            % Start Simulation
            bits = randi([0,1], 1, nSym*log2(Mary));
            
            if Mary == 2
                msg = bits;
            else
                msg = reshape(bits,[log2(Mary),nSym]);
                msg = (bi2de(msg.','left-msb'))';
            end
            
            msg1 = msg(1:length(msg)/2);
            msg2 = msg((length(msg)/2)+1:length(msg));
            
            s1 = qammod(msg1,Mary,'UnitAveragePower',true);
            s2 = qammod(msg2,Mary,'UnitAveragePower',true);
            x = [s1; s2];
            
            n1 = 10^(-EbNo(ebno_index)/20)*(1/sqrt(2)*(randn(1,length(s1)) + 1i*randn(1,length(s1))));
            n2 = 10^(-EbNo(ebno_index)/20)*(1/sqrt(2)*(randn(1,length(s2)) + 1i*randn(1,length(s2))));
            n =[n1; n2];
            
            % Send through channel
            % [r1 ; r2] = [h11 h12; h21 h22][s1; s2] + [n1; n2]
            y = H*x + n;
            
            % Zero-Forcing
            Wzf = inv(H);
            y_hat = Wzf * y;
            
            % Demodulate
            y1 = y_hat(1,:);
            y2 = y_hat(2,:);
            
            rx1 = qamdemod(y1,Mary,'UnitAveragePower',true);
            rx2 = qamdemod(y2,Mary,'UnitAveragePower',true);
            
            rx = [rx1 rx2];
            
            %Symbols back to bits
            if Mary == 2
                rxMSG = rx;
            else
                rxMSG = de2bi(rx.','left-msb')';
                rxMSG = reshape(rxMSG, [1,nSym*log2(Mary)]);
            end
            
            [~,ber(i,j,k)] = biterr(bits,rxMSG);
            
        end
    end
end
% Average trials
ber_table = mean(ber,3);

format short e
str=['BER table: ZF, EbNo(dB) = ', num2str(EbNo(ebno_index))]
T = array2table(ber_table,'RowNames',{'fm = 10Hz','fm = 100Hz','fm = 1kHz'},'VariableNames',{'BSPK','4 QAM','16 QAM','64 QAM'})

clear all;

%% Part 1b: 2x2 MIMO System with MMSE Channel Estimation

% This program simulates a 2x2 MIMO system with 1000 symbols. The signal is
% modulated with BPSK, 4 QAM, 16 QAM, and 64 QAM to compare BER and data
% rates. The signal is passed through three different flat fading Rayleigh 
% channels, generated with max doppler shifts of 10, 100 and 1000 Hz. The
% MMSE detection algorithm is implemented to minimize ISI and the effect of noise. 
% Assuming channel state information at the reciever (CSIR), a combination of the 
% inverse of the channel gain matrix and variance of the noise is applied to the recieved signal. 
% The simulation is repeated 10,000 times, and the trials are averaged. 
% My goal was to obtain a bit error rate on the order of 1e-3 
% for an SNR of 20 and doppler frequency shift of 10 Hz.
% 4 QAM reached the lowest BER, with a data rate of 2 bits per symbol.

% Parameters
N = 2^0; % number of frequency domain points for Rayleigh Channel generation
nSym = 1000; % number of symbols per packet
numIter = 10000; % number of iterations

EbNo = 0:25;
ebno_index = 21;
SNR = 20;

fm = [10 100 1000]; % max doppler frequency shift (in Hz)
M = [2 4 16 64];

for i = 1:length(fm)% LOOP 1 for Max Doppler Frequency Shift
    for j = 1:length(M)% LOOP 2 for M-ary number
        for k = 1:numIter% LOOP 3 for Number of Iterations

            % Get Rayleigh Channel Coefficient Matrix
            h11 = raychan(fm(i),N);
            h12 = raychan(fm(i),N);
            h21 = raychan(fm(i),N);
            h22 = raychan(fm(i),N);
            H = [h11 h12; h21 h22];
            
            % Set QAM
            Mary = M(j);
            
            % Start Simulation
            bits = randi([0,1], 1, nSym*log2(Mary));
            
            % Bits to symbols
            if Mary == 2
                msg = bits;
            else
                msg = reshape(bits,[log2(Mary),nSym]);
                msg = (bi2de(msg.','left-msb'))';
            end
            
            % 2 Signals
            msg1 = msg(1:length(msg)/2);
            msg2 = msg((length(msg)/2)+1:length(msg));
            
            % Modulate
            s1 = qammod(msg1,Mary,'UnitAveragePower',true);
            s2 = qammod(msg2,Mary,'UnitAveragePower',true);
            x = [s1; s2];
            
            % Create noise 
            n1 = 10^(-EbNo(ebno_index)/20)*(1/sqrt(2)*(randn(1,length(s1)) + 1i*randn(1,length(s1))));
            n2 = 10^(-EbNo(ebno_index)/20)*(1/sqrt(2)*(randn(1,length(s2)) + 1i*randn(1,length(s2))));
            n =[n1; n2];
            
            % Send signal through rayleigh channel and add noise
            % [r1 ; r2] = [h11 h12; h21 h22][s1; s2] + [n1; n2]
            y = H*x + n;
            
            % MMSE to eliminate ISI
            No = [var(n1) var(n2)];
            H_H = transpose(conj(H));% H^H
            
            Wmmse = inv(H_H*H + No*eye(2))*H_H;
            
            y_hat = Wmmse * y;
            
            % Demodulate
            y1 = y_hat(1,:);
            y2 = y_hat(2,:);
            
            rx1 = qamdemod(y1,Mary,'UnitAveragePower',true);
            rx2 = qamdemod(y2,Mary,'UnitAveragePower',true);
            
            rx = [rx1 rx2];
            
            %Symbols back to bits
            if Mary == 2
                rxMSG = rx;
            else
                rxMSG = de2bi(rx.','left-msb')';
                rxMSG = reshape(rxMSG, [1,nSym*log2(Mary)]);
            end
            
            [~,ber(i,j,k)] = biterr(bits,rxMSG);
            
        end% LOOP 3 for Number of Iterations
    end% LOOP 2 for M-ary number
end% LOOP 1 for Max Doppler Frequency Shift

% Average trials
ber_table = mean(ber,3);

format short e
title = ['BER table: MMSE, EbNo(dB) = ', num2str(EbNo(ebno_index))]
Table = array2table(ber_table,'RowNames',{'fm = 10Hz','fm = 100Hz','fm = 1kHz'},'VariableNames',{'BSPK','4 QAM','16 QAM','64 QAM'})

clear all;
%% Part 1c: 2x2 MIMO System with Pre-coding

% This program simulates a 2x2 MIMO system with 1000 symbols. The signal is
% modulated with BPSK, 4 QAM, 16 QAM, and 64 QAM to compare BER and data
% rates. The signal is passed through three different flat fading Rayleigh 
% channels, generated with max doppler shifts of 10, 100 and 1000 Hz. 
% Assuming perfect channel knowledge at the transmitter (CSIT), precoding is
% implemented to minimize ISI by following the method described in "Broadband
% MIMO-OFDM Wireless Communications" by Stuber et al. The orthogonal SVD
% decomposition of the channel gain matrix is obtained, using V to
% prefilter an U to postfilter.
% The simulation is repeated 10000 times, and the trials are averaged. 
% My goal was to obtain a bit error rate on the order of 1e-4 
% for an SNR of 20. 4 QAM reached the lowest BER, with a data rate of 2 bits per symbol.

% Parameters
N = 2^0; % number of frequency domain points for Rayleigh Channel generation
nSym = 1000; % number of symbols per packet
numIter = 10000; % number of iterations

EbNo = 0:25;
ebno_index = 21;
SNR = 20;

fm = [10 100 1000]; % max doppler frequency shift (in Hz)
M = [2 4 16 64];

for i = 1:length(fm)% LOOP 1 for Max Doppler Frequency Shift
    for j = 1:length(M)% LOOP 2 for M-ary number
        for k = 1:numIter% LOOP 3 for Number of Iterations

            % Get Rayleigh Channel Coefficient Matrix
            h11 = raychan(fm(i),N);
            h12 = raychan(fm(i),N);
            h21 = raychan(fm(i),N);
            h22 = raychan(fm(i),N);
            H = [h11 h12; h21 h22];
            
            % Get SVD decomposition
            [U,S,V] = svd(H);
            
            % Set M-ary for QAM
            Mary = M(j);
            
            % Start Simulation
            bits = randi([0,1], 1, nSym*log2(Mary));
            
            if Mary == 2
                msg = bits;
            else
                msg = reshape(bits,[log2(Mary),nSym]);
                msg = (bi2de(msg.','left-msb'))';
            end
            
            msg1 = msg(1:length(msg)/2);
            msg2 = msg((length(msg)/2)+1:length(msg));
            
            s1 = qammod(msg1,Mary,'UnitAveragePower',true);
            s2 = qammod(msg2,Mary,'UnitAveragePower',true);
            x_hat = [s1; s2];
            
            % Pre-filter
            x = V*x_hat;
            
            % Filter with Rayleigh Channel and add noise
            y = awgn(H*x,SNR+10*log10(log2(Mary)));
            
            % Post-filter
            y_hat = U'*y;
            y_hatnoise = inv(S)*y_hat; 
            
            y1 = y_hatnoise(1,:);
            y2 = y_hatnoise(2,:);
            
            % Demodulate
            y1 = y_hat(1,:);
            y2 = y_hat(2,:);
            
            rx1 = qamdemod(y1,Mary,'UnitAveragePower',true);
            rx2 = qamdemod(y2,Mary,'UnitAveragePower',true);
            
            rx = [rx1 rx2];
            
            %Symbols back to bits
            if Mary == 2
                rxMSG = rx;
            else
                rxMSG = de2bi(rx.','left-msb')';
                rxMSG = reshape(rxMSG, [1,nSym*log2(Mary)]);
            end
            
            [~,ber(i,j,k)] = biterr(bits,rxMSG);
            
        end
    end
end
% Average Trials 
ber_table = mean(ber,3);

format short e
str =['BER: precoding, SNR(dB) = ', num2str(SNR)]
T = array2table(ber_table,'RowNames',{'fm = 10Hz','fm = 100Hz','fm = 1kHz'},'VariableNames',{'BSPK','4 QAM','16 QAM','64 QAM'})


%% Functions
function [h] = raychan(fm, N)
% Purpose: Generate a Rayleigh Channel
% Inputs: fm == maximum Doppler frequency shift; N == number of frequency domain points

fc = 0; % @ baseband
f = linspace(-fm*0.9999,fm*0.9999,N); % f = fc +/- fm
% Fading Spectrum
SE = 1.5./(pi*fm*(sqrt(1-((f-fc)/fm).^2)));
%plot(f,SE) % Doppler PSD

delta_f = 2*fm/(N-1); % spacing between spectral lines
T = 1/delta_f;

% Gaussians Samples -- In-Phase
gaussI_pos = (1/sqrt(2))*(randn(1,length(N/2)) + 1i.*randn(1,length(N/2)));
gaussI_neg = fliplr(conj(gaussI_pos));
gaussI = [gaussI_neg gaussI_pos];

% Gaussian Samples -- Quadrature
gaussQ_pos = (1/sqrt(2))*(randn(1,length(N/2)) + 1i.*randn(1,length(N/2)));
gaussQ_neg = fliplr(conj(gaussQ_pos));
gaussQ = [gaussQ_neg gaussQ_pos];

% Multiply by fading spectrum
rayI_freq = gaussI.*sqrt(SE);
rayQ_freq = gaussQ.*sqrt(SE);

% Convert to time domain with IFFT
rayI_time = ifft(rayI_freq,N);
rayQ_time = ifft(rayQ_freq,N);

% N-point time series of a simulated Rayleigh fading signal
sumofsquares = rayI_time.^2 + rayQ_time.^2;
h = sqrt(sumofsquares);
end

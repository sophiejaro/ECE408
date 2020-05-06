clc
close all
clear all

% The purpose of this script is to create BER performance comparison of
% coherent BPSK with MRRC and two-branch transmit diversity in Rayleigh
% fading

% Plot with
%     -no diversity
%     -MRRC(1 Tx, 2 Rx)
%     -MRRC(1 Tx, 4 Rx)
%     -Alamouti Scheme(2 Tx, 1 Rx)
%     -Alamouti Scheme(2 Tx, 2 Rx)

EbNo = 0:2.5:50;
messageLength = 1e4;

% no diversity
M = 1; % 1 recieve antenna
BER_vec0 = MRRCscheme(M,EbNo,messageLength);
semilogy(EbNo,BER_vec0,'b:')
xlim([1 50])
ylim([1e-6 1])
xlabel('Eb/No (dB)');
ylabel('Pb, Bit Error Rate (BER)');

hold on

% MRRC(1 Tx, 2 Rx)
M = 2; % 2 recieve antenna
BER_vec1 = MRRCscheme(M,EbNo,messageLength);
semilogy(EbNo,BER_vec1,'r--')

hold on

% MRRC(1 Tx, 4 Rx)
M = 4; % 4 recieve antennas
BER_vec2 = MRRCscheme(M,EbNo,messageLength);
semilogy(EbNo,BER_vec2,'b--')

hold on

% Alamouti Scheme(2 Tx, 1 Rx)
M = 1; % 1 recieve antenna
BER_vec3 = alamoutiScheme(M,EbNo,messageLength);
semilogy(EbNo,BER_vec3,'r-')

hold on

% Alamouti Scheme(2 Tx, 2 Rx)
M = 2; % 2 recieve antenna
BER_vec4 = alamoutiScheme(M,EbNo,messageLength);
semilogy(EbNo,BER_vec4,'b-')

legend('no diversity scheme', 'MRRC (1Tx, 2Rx)', 'MRRC (1Tx, 4Rx)', 'Alamouti (2Tx, 1Rx)', 'Alamouti (2Tx, 2Rx)')

function [BER_vec] = MRRCscheme(M,EbNo,messageLength)

% Inputs
% M is number of recieve antennas (1 or 2)
% EbNo is the energy per bit over noise power vector

% Parameters
N = 2; % number of transmit antennas

for ebno_index = 1:length(EbNo)
    
    % Generate signal
    data = randi([0,1],messageLength,1);
    
    % BPSK modulate
    tx = pskmod(data,2);
    
    % Replicate for all recieve antennas
    tx_M = tx(:,ones(1,M));
    
    % Define Flat Fading Rayleigh Filter
    H =(1/sqrt(2))*(randn(messageLength, M) + 1i*randn(messageLength, M));
    tx_M_filt = H.*tx_M;
    
    % Add Noise
    rx = awgn(tx_M_filt, EbNo(ebno_index));
    
    % Combiner
    for i = 1:M
        s(:,i) = rx(:,i).* conj(H(:, i));
    end
    
    % Maximum Likelihood Detector
    for i = 1:length(s)
        x = s(i,:);
        y1 = complex(1);
        y2 = complex(-1);
        if (sum((x-y1).^2).^0.5) <= (sum((x-y2).^2).^0.5)
            si(i,:) = y1;
        else
            si(i,:) = y2;
        end
    end
    
    % BPSK demodulate
    out = pskdemod(si,2);
    
    % Calculate BER
    [~,b] = biterr(out, data);
    
    BER_vec(ebno_index,:) = b;
end
end

function [BER_vec] = alamoutiScheme(M,EbNo,messageLength)

% Inputs
% M is number of recieve antennas (1 or 2)
% EbNo is the energy per bit over noise power vector

% Parameters
N = 2; % number of transmit antennas

%Initialize arrays
txEnc = zeros(messageLength,N);
H = zeros(messageLength,N,M);
rx = zeros(messageLength,M);
s0_tilda = zeros(messageLength/N,M);
s1_tilda = s0_tilda;
s = zeros(messageLength,M);
si = zeros(size(s));

for ebno_index = 1:length(EbNo)
    
    % Generate signal
    data = randi([0,1],messageLength,1);
    
    % BPSK modulate
    tx = pskmod(data,2);
    
    % Alamouti Space-Time Encoding and Transmission Sequence
    % txEnc = [s0 s1; -s1* s0*]
    s0 = tx(1:2:end);
    s1 = tx(2:2:end);
    txEnc(1:2:end,:) = [s0 s1];
    txEnc(2:2:end,:) = [-conj(s1) conj(s0)];
    
    % Define Rayleigh Fading Channel h = [h0 h1]
    H(1:2:end,:,:) = (randn(messageLength/2,N,M) + (1/sqrt(2))*1i*randn(messageLength/2,N,M)); % defines h0
    H(2:2:end,:,:) = H(1:2:end,:,:); % defines h1
    
    % Recieved Signal
    % r0 = h0s0 + h1s1 + n0
    % r1 = -h0(s1*) + h1(s0*) + n1
    for i = 1:M
        rx(:,i) = awgn(sum(H(:,:,i).*txEnc,2)/sqrt(N), EbNo(ebno_index));
    end
    
    % Combining Scheme
    % s0_tilda = (h0*)r0 + h1(r1*);
    % s1_tilda = (h1*)r0 - h0(r1*)
    h_index = 1:2:length(H);
    for i = 1:M
        s0_tilda(:,i) = (conj(H(h_index,1,i))).*(rx(1:2:end,i)) + H(h_index,2,i).*conj(rx(2:2:end,i));
        s1_tilda(:,i) = (conj(H(h_index,2,i))).*(rx(1:2:end,i)) - H(h_index,1,i).*conj(rx(2:2:end,i));
    end
    s(1:2:end,:) = s0_tilda;
    s(2:2:end,:) = s1_tilda;
    
    % Maximum Likelihood Detector
    for i = 1:length(s)
        x = s(i,:);
        y1 = complex(1);
        y2 = complex(-1);
        if (sum((x-y1).^2).^0.5) <= (sum((x-y2).^2).^0.5)
            si(i,:) = y1;
        else
            si(i,:) = y2;
        end
    end
    %si = complex(si);
    
    % BPSK demodulate
    out = pskdemod(si,2);
    
    % Calculate BER
    [~,b] = biterr(out, data);
    
    BER_vec(ebno_index,:) = b;
end
end

function [H] = rayleigh(N,M)
fm = N/2;
fc = 0;

for i = 1:M
    gauss1 = randn(M,2*fm) + 1i.*randn(M,2*fm);
    gauss1(:,1:fm) = conj(gauss1(:,fm+1:end));
    
    gauss2 = randn(1,2*fm) + 1i.*randn(1,2*fm);
    gauss2(:,1:fm) = conj(gauss2(:,fm+1:end));
    
    f = linspace(-fm,fm,N);
    SE = 1.5./(pi.*fm.*sqrt(1-((f-fc)/fm).^2));
    
    m1 = gauss1.*SE;
    m2 = gauss2.*SE;
    
    m1(:,1)=0;
    m1(:,end)=0;
    
    m2(:,1)=0;
    m2(:,end)=0;
    
    m1m = ifft(m1);
    m2m = ifft(m2);
    
    H(i) = ((m1m.^2 + m2m.^2).^(1/2))';
end
end

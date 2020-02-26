% Wireless Communications
% Sophie Jaro
% 02/19/20

% LTE Standard Simulaion
% This project models the LTE 3gpp wireless standard downlink reciever and
% transmitter system. Particularly, it attempts to implement
% OFDMA (Orthogonal Frequency Division Multiple Access), the basis of the
% LTE PHY layer. The model simulates 1 frame of data (10 ms) consisting of 10
% subframes (1 ms each). Each subframe contains 2 slots, each with 14 symbols.
% Each resource block is modulated with QAM. OFDM processing of each symbol follows.
% Guard bands are added and the symbol is sent through the IFFT. Cyclic prefixes 
% are prepended. At this point, the peak power is computed and plotted. 
% The symbol is then passed through a noisy channel. To recover the signal,
% the cyclic prefixes are removed, then the symbol is passed through the FFT,
% and the guard bands are removed. The data block is then demodulated and a
% BER is computed.

clear all;
close all;
clc

% Parameters %

% Bandwidth = 5 MHz 
    % LTE also supports channel bandwidth of 1.4MHz, 3MHz, 5MHz, 10MHz, 15MHz, and 20MHz 
% Resource Blocks = 25
% Number of subcarriers = 512
nFFT = 512;
nSC = 600;

% Cyclic Prefix: Normal

% This program can simulate either 16 or 64 QAM. Other modulation schemes
% supported by LTE are QPSK, QAM, and higher orders of QAM
modulation = input('Choose [16] or [64] QAM? \n');
if modulation == 16
    Mary = 16;
    bps = 4; % 4 bits per symbol for 16 QAM
    modulationString = '16 QAM';
elseif modulation == 64
    Mary = 64;
    bps = 6; % 4 bits per symbol for 64 QAM
    modulationString = '64 QAM';
end
    
SNR = input('Choose SNR Value: [0] or [10] \n');

% Definitions %
% QPSK Table
qpsk_table = 1/sqrt(2)*[1+j, 1-j,-1+j,-1-j];

% Average Power
avgPower = 1.3778e-004;

% LTE Frame Structure Type 1: FDD (Frequency Division Duplexing)
% one frame is 10ms, one subframe is 1ms, one slot is 0.5 ms
frames = 1;
spf = 10; % subframes per frame
sps = 2; % slots per subframe
OFDMps = 7; % OFDM symbols per slot

nTTI = 10; % number of Transmission Time Intervals (1ms), 10 TTI = 1 frame
count = 1;

%% Begin Simulation 
for nSubframes = 1:nTTI
    txDataSym = [];
    bitsMatrix = [];
    for nSubcarriers = 1:nSC
        bits = randi([0,1],sps*OFDMps*bps,1);
        bitsMatrix(:,nSubcarriers) = bits;
        nSym = length(bits)/bps;
        data = reshape(bits,[log2(Mary),nSym]);
        data = (bi2de(data.','left-msb'))';
        dataModulated = qammod(data,Mary);
        txDataSym = [txDataSym; dataModulated];
    end
    
    rxDataSym=[];
    % There are 14 data symbols per subframe (2 slots with 7 OFDM symbols
    % per slot)
    for symbol = 1:14
        % Determine cyclic prefix length by symbol 
        % If pilots bits were added, it would be in this step
        if (symbol == 1 || symbol == 8) % Add long CP (160)
            CP = 159;
            txSymbol = txDataSym(:,symbol); % Each column represents a symbol
%             txSymbol(1:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1)); % pilot
%             txSymbol(7:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1)); 
        elseif (symbol == 5 || symbol == 12) % Add normal CP (144), insert pilot
            txSymbol = txDataSym(:,symbol); % Each column represents a symbol
%             txSymbol(4:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1));
%             txSymbol(10:12:end) = sqrt(2)*qpsk_table(randi([1,4], 50, 1));
            CP = 143; % Add normal CP (144)
        else
            txSymbol = txDataSym(:,symbol); % Each column represents a symbol
            CP = 143; % Add normal CP (144)
        end
        
        % Add guards, center around 0 DC
        txGuards = [zeros(723,1); txSymbol(1:300); 0;...
            txSymbol(301:end); zeros(724,1)];
            % The txGuard is plotted as a scatterplot below
            
        % Apply IFFT
        txOFDM = ifft(txGuards); % Frequency to time domain
             % The signal in the time domain is plotted below
             
        % Insert cyclic prefix at beginning of symbol
        % each OFDM symbol is preceded by a cyclic prefix to eliminate ISI
        txOFDMcp = txOFDM;
        txOFDMcp = [txOFDMcp(end-CP:end); txOFDMcp];
        
            % Compute peak power
            peakPower = max(abs(txOFDMcp).^2); 
            PAPR(count) = peakPower/avgPower;
            count = count + 1;
        
        % AWGN
        txNoisy = awgn(txOFDMcp,SNR+10*log10(log2(Mary)), 'measured');

        % Remove cyclic prefix from symbol
        rxRCP = txNoisy(CP+2:end);
        
        % Apply FFT
        rxFFT = fft(rxRCP); % Time to frequency domain
        
        % Remove guards
        rxRG = rxFFT(724:723+300); 
        rxRG2 = rxFFT(724+301:724+600);
        rxSymbol = [rxRG; rxRG2];
        rxDataSym(:,symbol) = rxSymbol;     
    end % symbol loop
     
   % Demodulation Loop
   rxBits = [];
   rxMSGMatrix = [];
   for nSubcarrier = 1:nSC
       rx = qamdemod(rxDataSym(nSubcarrier,:),Mary);
       rxMSG = de2bi(rx.','left-msb')';
       rxMSG = reshape(rxMSG, [1,nSym*log2(Mary)])';
       rxMSGMatrix(:,nSubcarrier) = rxMSG;
   end
   [~, berVec(nSubframes)] = biterr(bitsMatrix, rxMSGMatrix);
    
end % end subframe loop

ber = mean(berVec',1)

%% Checkpoints
scatterplot(txGuards)
titleString = sprintf('Input to IFFT, txGuards: %s Modulated with DC at 0',...
    modulationString);
title(titleString)

scatterplot(txSymbol)
titleString = sprintf('txSymbol: %s Modulated without DC', modulationString);
title(titleString)

scatterplot(txOFDM)
title('Output of IFFT')

scatterplot(rxFFT)
title('Recieved signal: Output of FFT')


%% Plot Complementary Cumulative Distribution Function of the Peak to Average Power Ratio
bins = sort(PAPR);
cdf = linspace(0,1,length(bins));
ccdf = 1 - cdf;

figure
ylim([1e-2 1])
semilogy(10*log10(bins), ccdf)
grid
xlabel('PAPR (dB)')
ylabel('Probability')
titleString = sprintf('PAPR CCDF for %s',modulationString); 
title(titleString);
legend(modulationString)

pause

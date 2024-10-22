% Skeleton code for simulation chain

% History:
%   2000-06-28  written /Stefan Parkvall
%   2001-10-22  modified /George Jongren

clear
% Initialization
EbN0_db = 0:10;                     % Eb/N0 values to simulate (in dB)
nr_bits_per_symbol = 2;             % Corresponds to k in the report
nr_guard_bits = 10;                 % Size of guard sequence (in nr bits)
                                    % Guard bits are appended to transmitted bits so
                                    % that the transients in the beginning and end
                                    % of received sequence do not affect the samples
                                    % which contain the training and data symbols.
nr_data_bits = 1000;                % Size of each data sequence (in nr bits)
nr_training_bits = 12;             % Size of training sequence (in nr bits)
nr_blocks = 50;                     % The number of blocks to simulate
Q = 8;                              % Number of samples per symbol in baseband

% Define the pulse-shape used in the transmitter. 
% Pick one of the pulse shapes below or experiemnt
% with a pulse of your own.
pulse_shape = ones(1, Q);
% pulse_shape = root_raised_cosine(Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.
nr_errors = zeros(1,length(EbN0_db));   % Error counter
nr_errors2 = zeros(1,length(EbN0_db));
phihat_r1 = zeros(1,length(EbN0_db));
phihat_r2 = zeros(1,length(EbN0_db));
t_psamp = zeros(1,length(EbN0_db));
t_psamp2 = zeros(1,length(EbN0_db));
% Loop over several blocks to get sufficient statistics.
for snr_point = 1:length(EbN0_db)
    for blk = 1:nr_blocks

    %%%
    %%% Transmitter
    %%%

    % Generate training sequence.
    b_train = training_sequence(nr_training_bits);
    
    % Generate random source data {0, 1}.
    b_data = random_data(nr_data_bits);

    % Generate guard sequence.
    b_guard = random_data(nr_guard_bits);
 
    % Multiplex training and data into one sequence.
    b = [b_guard b_train b_data b_guard];
    
    % Map bits into complex-valued QPSK symbols.
    d = qpsk(b);

    % Upsample the signal, apply pulse shaping.
    tx = upfirdn(d, pulse_shape, Q, 1);
    
    %%% AWGN Channel
    % Compute variance of complex noise according to report.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);
    
    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(tx))+1j*randn(size(tx)));

    % Received signal.
    rx = tx + n;

    %%% Receiver  
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
    
    % Synchronization. The position and size of the search window
    % is here set arbitrarily. Note that you might need to change these
    % parameters. Use sensible values (hint: plot the correlation
    % function used for syncing)! 
    t_start = 1+Q*nr_guard_bits/2;
    t_end = t_start+50;
    t_samp = sync(mf, b_train, Q, t_start, t_end);
    t_psamp(snr_point) = t_samp - 48;
    
    % Down sampling. t_samp is the first sample, the remaining samples are all
    % separated by a factor of Q. Only training+data samples are kept.
    r = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);

    % Phase estimation and correction.
    phihat = phase_estimation(r, b_train);
    phihat_r1(snr_point) = phihat;
    r = r * exp(-1i*phihat);
        
    % Make decisions. Note that dhat will include training sequence bits
    % as well.
    bhat = detect(r);
    
    % Count errors. Note that only the data bits and not the training bits
    % are included in the comparison. The last data bits are missing as well
    % since the whole impulse response due to the last symbol is not
    % included in the simulation program above.
    temp=bhat(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors(snr_point) = nr_errors(snr_point) + sum(temp);
    
    %%% two-path ISI Channel
    tau = 8;
    tx2 = multipath(tx,tau);
    %%% AWGN Channel
    % Compute variance of complex noise according to report.
    sigma_sqr2 = norm(pulse_shape)^2 / nr_bits_per_symbol / 10^(EbN0_db(snr_point)/10);
    % Create noise vector.
    n2 = sqrt(sigma_sqr2/2)*(randn(size(tx2))+1j*randn(size(tx2)));
    rx2 = tx2 + n2;
    % Matched filtering.
    mf2=conv(mf_pulse_shape,rx2);
    t_samp2 = sync(mf2, b_train, Q, t_start, t_end);
    t_psamp2(snr_point) = t_samp2 - 48;
    r2 = mf2(t_samp2:Q:t_samp2+Q*(nr_training_bits+nr_data_bits)/2-1);
    phihat2 = phase_estimation(r2, b_train);
    phihat_r2(snr_point) = phihat2;
    r2 = r2 * exp(-1i*phihat2);
    bhat2 = detect(r2);
    temp2 = bhat2(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors2(snr_point) = nr_errors2(snr_point) + sum(temp2);
    % Next block.
    end
end

% Compute the BER. 
BER_AWGN = nr_errors / nr_data_bits / nr_blocks;
BER_ISI = nr_errors2 / nr_data_bits / nr_blocks;
figure(1);
plot(EbN0_db,BER_AWGN);
hold on;
plot(EbN0_db,BER_ISI,'r');
set(gca, 'YScale', 'log');
legend('AWGN channel','two-path ISI channel');
xlabel('SNR/dB');
ylabel('BER');

figure(2);
plot(EbN0_db,phihat_r1);
hold on
plot(EbN0_db,phihat_r2,'r');
legend('AWGN channel','two-path ISI channel');

figure(3);
plot(EbN0_db,t_psamp);
hold on
plot(EbN0_db,t_psamp2,'r');
legend('AWGN channel','two-path ISI channel');

% plot eyediagram
 figure(4);
eyediagram(r,4);
 figure(5);
 eyediagram(r2,4);
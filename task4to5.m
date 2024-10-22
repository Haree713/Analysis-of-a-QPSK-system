% Skeleton code for simulation chain
% Simulation for task 4 and 5 (PSD and eye diagram)

clear

% Initialization
EbN0_db = 0:10;                    
nr_bits_per_symbol = 2;            
nr_guard_bits = 10;                                                                                        
nr_data_bits = 1000;              
nr_training_bits = 100;             
nr_blocks = 50;                   
Q = 8;                              

% Define the pulse-shape used in the transmitter. 

  pulse_shape = ones(1, Q);
% pulse_shape = root_raised_cosine(Q);

% Matched filter impulse response. 
mf_pulse_shape = fliplr(pulse_shape);


% Loop over different values of Eb/No.

nr_errors1 = zeros(1,length(EbN0_db));  
nr_errors2 = zeros(1,length(EbN0_db));
phihat_r1 = zeros(1,length(EbN0_db));
phihat_isi = zeros(1,length(EbN0_db));
t_psamp1 = zeros(1,length(EbN0_db));
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

    %%%
    %%% two path isi
    %%%

    t = 6;
    tx2 = multipath(tx,t);

    % Compute variance of complex noise
    % The code is two separate lines due to needed space.
    sigma_sqr2 = norm(pulse_shape)^2 / nr_bits_per_symbol /10^(EbN0_db(snr_point)/10);

    % Create noise vector.
    n2 = sqrt(sigma_sqr2/2)*(randn(size(tx2))+1j*randn(size(tx2)));
    rx2 = tx2 + n2;


    %%%
    %%% AWGN Channel
    %%%
    
    % Compute variance of complex noise
    % The code is two separate lines due to needed space.
    sigma_sqr = norm(pulse_shape)^2 / nr_bits_per_symbol /10^(EbN0_db(snr_point)/10);

    % Create noise vector.
    n = sqrt(sigma_sqr/2)*(randn(size(tx))+j*randn(size(tx)));

    % Received signal.
    rx = tx + n;

    %%%
    %%% Receiver
    %%%
    
    % Matched filtering.
    mf=conv(mf_pulse_shape,rx);
    mf2=conv(mf_pulse_shape,rx2);
    
    % Synchronization
  
    t_start=1+Q*nr_guard_bits/2;
    t_end=t_start+50;
    t_samp = sync(mf, b_train, Q, t_start, t_end);
    t_samp2 = sync(mf2, b_train, Q, t_start, t_end);
    t_psamp2(snr_point) = t_samp2 - 48;
    
    % Down sampling

    r = mf(t_samp:Q:t_samp+Q*(nr_training_bits+nr_data_bits)/2-1);
    r2 = mf2(t_samp2:Q:t_samp2+Q*(nr_training_bits+nr_data_bits)/2-1);

    % Phase estimation and correction

    phihat = phase_estimation(r, b_train);
    r = r * exp(-j*phihat);

     phihat2 = phase_estimation(r2, b_train);
    phihat_isi(snr_point) = phihat2;
     r2 = r2 * exp(-1i*phihat2);
        
    % decisions
    bhat = detect(r);
    bhat2 = detect(r2);
    
    % Count errors

    temp=bhat(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors1(snr_point) = nr_errors1(snr_point) + sum(temp);

    temp2 = bhat2(1+nr_training_bits:nr_training_bits+nr_data_bits) ~= b_data;
    nr_errors2(snr_point) = nr_errors2(snr_point) + sum(temp2);

    % Next block.
  end

  % Next Eb/No value.
end

% Compute the BER. 
BER_AWGN = nr_errors1 / nr_data_bits / nr_blocks;
BER_ISI = nr_errors2 / nr_data_bits / nr_blocks;
EbN0 = 10.^(EbN0_db/10);

figure(1);

plot(EbN0_db,BER_ISI,'r');
set(gca, 'YScale', 'log');
xlabel('SNR/dB');
ylabel('BER');
title('BER of Two-path ISI channel');

figure(2);

plot(EbN0_db,phihat_isi,'r');
title('Phase Accuracy of Two-path ISI channel');
xlabel('SNR/db');
ylabel('Phase error')

figure(3);

plot(EbN0_db,t_psamp2,'r');
title('Synchronization Accuracy of Two-path ISI channel for training bits = 8');
xlabel('SNR/db');
ylabel('Synchronization error');
 
figure;
periodogram(tx);

figure;
eyediagram(rx,4);    
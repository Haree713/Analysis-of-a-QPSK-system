function phihat = phase_estimation(r, b_train)

phihat = 0;
qpsk_train = qpsk(b_train);
r_train = r(1:length(qpsk_train));
min = norm(r_train - qpsk_train); %minimizing the norm
for phase = -pi:0.01:pi
    r_phi = r_train * exp(-1j*phase);
    min_length = norm(r_phi - qpsk_train);
    if min_length < min
        min = min_length;
        phihat = phase;
    end
end
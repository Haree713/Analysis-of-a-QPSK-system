function t_samp = sync(mf, b_train, Q, t_start, t_end)

% function t_samp = sync(mf, b_train, Q, t_start, t_end)
    qpsk_train = qpsk(b_train);
    len = length(qpsk_train);

%     Initialize the cross-correlation results array
    r = zeros(1, t_end - t_start + 1);

%     Calculate cross-correlation for each time shift
    for shift = t_start:t_end
%         Extract the matched filter output for the corresponding time shift
        a = mf(shift + (0:len-1) * Q);
        
%         Calculate the absolute value of the inner product
        r(shift - t_start + 1) = abs(sum(conj(a) .* qpsk_train));
    end

%     Find the index of the maximum cross-correlation value
    [~, t_samp] = max(r);
    t_samp = t_samp + t_start - 1;
end
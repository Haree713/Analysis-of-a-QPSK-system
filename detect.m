 function bhat = detect(r)
% bhat = detect(r)

n = length (r);
bhat = zeros(1,2*n);
for i=1:n
    if real(r(i)) > 0
        bhat(2*i-1) = 0;
    else
        bhat(2*i-1) = 1;
    end
    if imag(r(i)) > 0;
        bhat(2*i) = 0;
    else
        bhat(2*i)= 1;
    end
end
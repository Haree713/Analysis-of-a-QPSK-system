function d = qpsk(b)

d=zeros(1,length(b)/2);
%definition of the QPSK symbols using Gray coding.
for n=1:length(b)/2
    p=b(2*n);
    imp=b(2*n-1);
    if (imp==0)&(p==0)
        d(n)=exp(j*pi/4);%45 degrees
    end
    if (imp==1)&(p==0)
        d(n)=exp(j*3*pi/4);%135 degrees
    end
    if (imp==1)&(p==1)
        d(n)=exp(j*5*pi/4);%225 degrees
    end
    if (imp==0)&(p==1)
        d(n)=exp(j*7*pi/4);%315 degrees
    end
end
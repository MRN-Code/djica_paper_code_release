function D = divisor(N)
K = 1:N;
D = K(rem(N,K)==0);
D = D(1:(end-1));
end

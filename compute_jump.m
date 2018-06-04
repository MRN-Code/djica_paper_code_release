function jump = brad_compute_jump(s,mmax,mmin,scale)
    if ~exist('mmax','var')
        mmax = 1;
    end
    if ~exist('mmax','var')
        mmin = Inf;
    end
    if ~exist('scale','var')
        scale = 1;
    end
    jump = min(max(floor(max(divisor(s))/scale),mmin),mmax);
    if isprime(s) && s > 2
        jump = min(max(floor(max(divisor(s-1))/scale),mmin),mmax);
    end
end
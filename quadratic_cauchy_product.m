function [ab] = quadratic_cauchy_product(a,b)
% Computes the Cauchy product of two sequences with possibly different lengths

m1 = length(a);
m2 = length(b);

ab = intval(0)*zeros(m1+m2-1,1);

for n=(1:(m1+m2-1))
    if (n < m1) 
        if (n < m2)
            ab(n) = sum(a(1:n).*b(n:-1:1));
        else
            ab(n) = sum(a((n-m2+1):n).*b(m2:-1:1));
        end
    else
        if (n < m2)
            ab(n) = sum(a(1:m1).*b(n:-1:(n-m1+1)));
        else
            ab(n) = sum(a((n-m2+1):m1).*b(m2:-1:(n-m1+1)));
        end
    end
end


end


function nphi = jingxuan(alpha)

fun = @(x) 1 ./ (1 + (x + 0.85).^30);

if alpha == 0;
    nphi = 1;
else
    nphi = integral(fun, 0, alpha) / alpha;
    
end
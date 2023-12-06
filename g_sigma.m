function result = g_sigma(sigma)
    result = sqrt(2*pi*sigma.^2) * exp((1+1/(4*sigma.^2)).^2 - 1);
end
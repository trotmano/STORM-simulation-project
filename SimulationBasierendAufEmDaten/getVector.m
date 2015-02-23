function vec = getVector(aoa, length, alpha)
    z = sin(aoa)*sin(alpha);
    y = sin(aoa)*cos(alpha);
    x = sqrt(1-y.^2 - z.^2);
    vec = [real(x),real(y),real(z)]*length;
end
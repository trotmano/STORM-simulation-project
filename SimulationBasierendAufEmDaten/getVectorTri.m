function vec = getVectorTri(aoa, length)
    alpha = rand(1)*2*pi;
    x = cos(aoa)*cos(alpha);
    z = sin(aoa);
    y = cos(aoa)*sin(alpha);
    
    vec = [real(x),real(y),real(z)]*length;
end
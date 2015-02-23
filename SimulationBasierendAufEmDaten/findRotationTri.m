function rotVec = findRotationTri(vec,normVec)
    unityVec = [0,1,0];
    targetVec = normVec;
    targetVec = targetVec / norm(targetVec);
    v= cross(unityVec,targetVec);
    s = norm(v);
    c = unityVec*targetVec';
    vx = [0,-v(3), v(2);v(3),0,-v(1);-v(2),v(1),0];
    R = [1,0,0;0,1,0;0,0,1]+vx+vx*vx*(1-c)/s^2;
    rotVec = R*vec'; %http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
%     
%     alpha = atan2(targetVec(2),targetVec(1));
%     drehy = [cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];
%     tmpV1 = drehy * vec';
%     beta = -acos(tmpV1' * targetVec'/(norm(tmpV1)*norm(targetVec)));
%     drehx=[1,0,0;0,cos(beta),-sin(beta);0,sin(beta),cos(beta)];
%     rotVec = drehx * tmpV1;
end

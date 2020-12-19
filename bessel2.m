a = 6378137; %WGS-84
b = 6356752.3142;
e2 = 0.00669437999013;
eps2 = 0.00673949674227;

%input B1 L1 B2 L2 output A1 A2 S
B1 = deg2rad(dms2degrees(B1));
L1 = deg2rad(dms2degrees(L1));
B2 = deg2rad(dms2degrees(B2));
L2 = deg2rad(dms2degrees(L2));

%1
U1 = atan(sqrt(1-e2)*tan(B1));
U2 = atan(sqrt(1-e2)*tan(B2));
l = L2-L1;
a1 = sin(U1)*sin(U2); a2 = cos(U1)*cos(U2);
b1 = cos(U1)*sin(U2); b2 = sin(U1)*cos(U2); %1.1
lambda = l;
Dlambda = 1;
while(Dlambda>0.001)
    lambda0 = lambda;
    p = cos(U2)*sin(lambda0);
    q = b1 - b2*cos(lambda0);
    A1 = abs(atan(p/q));
    if(p>0 && q>0) %A1=A1
    end
    if(p>0 && q<0)
        A1 = pi-A1;
    end
    if(p<0 && q<0)
        A1 = pi+A1;
    end
    if(p<0 && q>0)
        A1 = 2*pi-A1;
    end
    sinsig = p*sin(A1)+q*cos(A1); cossig = a1+a2*cos(lambda0);
    sigma = abs(atan(sinsig/cossig));
    if(cossig>0) %sigma=sigma
    end
    if(cossig<0)
        sigma = pi-sigma;
    end
    sinA0 = cos(U1)*sin(A1);
    sigma1 = atan(tan(U1)*sec(A1));
    cosA0 = sqrt(1-sinA0^2);
    e4 = e2^2; e6 = e2^3;
    kp2 = e2*cosA0^2; kp4 = kp2^2;
    alpha1 = (e2/2+e4/8+e6/16) - e2*(1+e2)*kp2/16 + 3*kp4*e2/128;
    beta1 = e2*(1+e2)*kp2/16 - e2*kp4/32;
    gamma1 = e2*kp4/256;
    xx = alpha1*sigma + beta1*sin(sigma)*cos(2*sigma1+sigma);
    xx = xx + gamma1*sin(2*sigma)*cos(4*sigma1+2*sigma);
    lambda = l+sinA0*xx;
    Dlambda = abs(lambda - lambda0)*206265;
end %1.2
%2
k2 = eps2*cosA0^2; k4 = k2^2; k6 = k2^3;
alpha = (1-k2/4 + 7*k4/64 - 15*k6/256)/b;
beta = k2/4 - k4/8 +37*k6/512;
gamma = k4/128-k6/128;
S = gamma*sin(2*sigma)*cos(4*sigma1+2*sigma);
S = (sigma-beta*sin(sigma)*cos(2*sigma1+sigma) - S) / alpha;
sinA2 = cos(U1)*sin(A1);
cosA2 = cos(U1)*cos(sigma)*cos(A1)-sin(U1)*sin(sigma);
tanA2 = sinA2/cosA2;
A2 = abs(atan(tanA2));
sinA1 = sin(A1);
if(sinA1<0 && tanA2>0) %A2=A2
end
if(sinA1<0 && tanA2<0)
    A2 = pi-A2;
end
if(sinA1>0 && tanA2>0)
    A2 = pi+A2;
end
if(sinA1>0 && tanA2<0)
    A2 = 2*pi-A2;
end

disp(rad2deg(A1));
disp(rad2deg(A2));
disp(S);

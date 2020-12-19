a = 6378137; %WGS-84
b = 6356752.3142;
e2 = 0.00669437999013;
eps2 = 0.00673949674227;

%input B1 L1 A1 S output B2 L2 A2
B1 = deg2rad(dms2degrees(B1));
L1 = deg2rad(dms2degrees(L1));
A1 = deg2rad(dms2degrees(A1));

%1
U1 = atan(sqrt(1-e2)*tan(B1)); %1.1
sinA0 = cos(U1)*sin(A1);
cosA0 = sqrt(1-sinA0^2);
sigma1 = atan(tan(U1)*sec(A1)); %1.2
k2 = eps2*cosA0^2; k4 = k2^2; k6 = k2^3;
alpha = (1-k2/4 + 7*k4/64 - 15*k6/256)/b;
beta = k2/4 - k4/8 +37*k6/512;
gamma = k4/128-k6/128;
sigma = alpha*S;
Dsigma = 1;
while(Dsigma>0.001)
    sigma0 = sigma;
    sigma = alpha*S+beta*sin(sigma0)*cos(2*sigma1+sigma0);
    sigma = sigma + gamma*sin(2*sigma0)*cos(4*sigma1+2*sigma0);
    Dsigma = abs(sigma - sigma0)*206265;
end %1.3
%2
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
end %2.1
sinU2 = sin(U1)*cos(sigma)+cos(U1)*cos(A1)*sin(sigma); %2.2
sinl = sin(A1)*sin(sigma); cosl = cos(U1)*cos(sigma)-sin(U1)*sin(sigma)*cos(A1);
tanl = sinl/cosl;
lambda = abs(atan(tanl));
if(tanl>0 && sinA1>0) %lambda=lambda
end
if(tanl<0 && sinA1>0)
    lambda = pi-lambda;
end
if(tanl<0 && sinA1<0)
    lambda = -lambda;
end
if(tanl>0 && sinA1<0)
    lambda = lambda-pi;
end %2.3
%3
B2 = atan(sinU2/sqrt(1-e2)/sqrt(1-sinU2^2)); %3.1
e4 = e2^2; e6 = e2^3;
kp2 = e2*cosA0^2; kp4 = kp2^2;
alpha1 = (e2/2+e4/8+e6/16) - e2*(1+e2)*kp2/16 + 3*kp4*e2/128;
beta1 = e2*(1+e2)*kp2/16 - e2*kp4/32;
gamma1 = e2*kp4/256;
xx = alpha1*sigma + beta1*sin(sigma)*cos(2*sigma1+sigma);
xx = xx + gamma1*sin(2*sigma)*cos(4*sigma1+2*sigma);
l = lambda-sinA0*xx;
L2 = L1+l; %3.2

disp(rad2deg(B2));
disp(rad2deg(L2));
disp(rad2deg(A2));

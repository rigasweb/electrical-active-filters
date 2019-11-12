%% Tzikas Tryfon Rigas
%% 8589

%% LowPass Inverse Chebyshev


clear

%initializations
a_1 = 8;
a_2 = 5;
a_3 = 8;
a_4 = 9;
m = 3;

f_p = 1.1 * ( 3 + m ) * 10^3 ;
f_s = 1.9 * f_p ;
a_min = 25 + ( max(1, a_3) - 5) * 3 / 4 ;
a_max = 0.55 + ( max(1,a_4) - 5) / 16 ;

w_p = 2 * pi * f_p;
w_s = 2 * pi * f_s;

%normalize the frequencies so that we have W_s=1
W_s = 1;
W_p = w_p / w_s;
n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(1/W_p));

n=ceil(n);

e =1 / sqrt(10^(a_min/10)-1);
a = (1/n)*asinh(1/e);
w_hp = 1/cosh((1/n)*acosh(1/e));
f_hp = w_hp * f_p;

%butterworth angles
y_k1 = 22.5;
y_k2 = -22.5;
y_k3 = 67.5;
y_k4 = -67.5;

%poles
p_k1 = -sinh(a)*cosd(y_k1) + 1i*cosh(a)*sind(y_k1);
p_k2 = real(p_k1) - 1i*imag(p_k1);
p_k3 = -sinh(a)*cosd(y_k3) + 1i*cosh(a)*sind(y_k3);
p_k4 = real(p_k3) - 1i*imag(p_k3);

%W_0 and Q for the poles
W_012 = (sqrt(real(p_k1)^2 + imag(p_k1)^2));
W_034 = (sqrt(real(p_k3)^2 + imag(p_k3)^2));

Q_12 = (sqrt(real(p_k1)^2 + imag(p_k1)^2))/abs(2*real(p_k1));
Q_34 = (sqrt(real(p_k3)^2 + imag(p_k3)^2))/abs(2*real(p_k3));

%inverse poles
InvW1 = 1 / W_012;
InvW2 = 1 / W_034;
S12 = -InvW1 / ( 2 * Q_12 );
S34 = -InvW2 / ( 2 * Q_34);
W12 = sqrt( InvW1^2 - S12^2 );
W34 = sqrt( InvW2^2 - S34^2 );
w_012 = InvW1;
w_034 = InvW2;


%ta mhdenika ths sunarthshs
w_z1 = sec( pi / (2 * n));
w_z3 = sec(3* pi / (2 * n));

w_0 = w_p /(((10^(a_max/10)-1))^(1/(2*n)));
f_0 = w_0/(2*pi);

%ta pragmatika metra twn polwn
w_12 = w_012 *w_s;
w_34 = w_034 *w_s;
%%
% a3=8 ara xrhsimopoioume kanoniko Notch kuklwma
% a4=9 ara theloume puknwth 1.0 ìF
% a2=5 ara kerdos 10 dB
%% Prwth Monada
W_z1 = w_z1 / w_034;  % W_z1>1
R11 = 1;
R12 = 4 * Q_34^2;
R13 = (W_z1^2)/(2*Q_34^2);
R14 = 1;
R15 = (4*Q_34^2)/(W_z1^2-1);
C1 = 1 / ( 2 * Q_34 ) ;
k_h1 = 1 / (1 + R13); %kerdos stis ypshles syxnothtes
k_l1 = k_h1 * (w_z1/w_034)^2; %keros stis xamhles syxnothtes

k_f1 = w_s * w_034;
k_m1 = 10^6 * C1 / k_f1;

%klimakopoihsh
R11 = R11 * k_m1;
R12 = R12 * k_m1;
R13 = R13 * k_m1;
R14 = R14 * k_m1;
R15 = R15 * k_m1;
C1 = 10^(-6);

%% Deuterh Monada
W_z3 = w_z3 / w_012;  % W_z3>1
R21 = 1;
R22 = 4 * Q_12^2;
R23 = (W_z3^2)/(2*Q_12^2);
R24 = 1;
R25 = (4*Q_12^2)/(W_z3^2-1);
C2 = 1 / ( 2 * Q_12 ) ;
k_h2 = 1 / (1 + R23); %kerdos stis ypshles syxnothtes
k_l2 = k_h2 * (w_z3/w_012)^2; %keros stis xamhles syxnothtes

k_f2 = w_s * w_012;
k_m2 = 10^6 * C2 / k_f2;

%klimakopoihsh
R21 = R21 * k_m2;
R22 = R22 * k_m2;
R23 = R23 * k_m2;
R24 = R24 * k_m2;
R25 = R25 * k_m2;
C2 = 10^(-6);


%% Transfer Function
T1 = tf([k_l2 0 k_l2*((w_z3*w_s)^2) ], [1, w_12/Q_12, w_12^2]);

T2 = tf([k_l1 0 k_l1*((w_z1*w_s)^2) ], [1, w_34/Q_34, w_34^2]);

T_all = T1 * T2;

%% Total Gain
k_total = k_l1 * k_l2;
%klimakopoihsh
gain = abs(evalfr(T_all, w_p * 1i));
k_desired = 10^0.5;
k = k_desired / gain;
%k<1
%ZA = R11 / k;
%ZB = R11 / (1-k); 
ZA = 1000;
ZB = 1000*(k_desired/k_total); 

T_total = k*T_all

plot_transfer_function(T1, [f_p f_s ])

plot_transfer_function(T2, [f_p f_s ])

plot_transfer_function(T_all, [f_p f_s ])

plot_transfer_function(T_total, [f_p f_s ])

ltiview({'bodemag'}, T1, T2, T_total)
 
InvSys_new = inv(T_all)
plot_transfer_function(InvSys_new, [f_p f_s ])

InvSys_new1 = inv(T_total)
plot_transfer_function(InvSys_new1, [f_p f_s ])

%% Fourier Analysis
fs= 2*10^3;
T=0.2;
dt=1/fs;
t=0:dt:(T-1/fs);
x = sawtooth(2*pi*20*t);

figure %plot tou simatos eisodou
plot(t,x)
title('Sawtooth wave') 
xlabel('t (sec)')
ylabel('Amplitude') 

%dimiourgia twn fasmatwn fourier: eisodou kai eksodou
N=T/dt;
xt=lsim(T_all,x,t);
figure
plot(t,xt)
n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=2*10^3*(0:(n/2))/n;
figure
plot(f,p1)
nfft=n;
y=fft(x,nfft);
figure
plot(xt)
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*fs/nfft; 
figure
plot(f,y_mag)

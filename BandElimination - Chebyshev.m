%% Tzikas Tryfon Rigas

%% BandElimination Chebyshev


clear

%initializations
a_1 = 8;
a_2 = 5;
a_3 = 8;
a_4 = 9;

f_0 = 2500;
f_1 = 1700 + 50 * a_3;
f_2 = f_0^2 / f_1;
D = (1/2.1) * (f_0^2-f_1^2) / f_1;
f_3 = (-D+sqrt(D^2+4*f_0^2))/2;
f_4 = f_0^2/f_3;

w_0 = 2 * pi * f_0;
w_1 = 2 * pi * f_1;
w_2 = 2 * pi * f_2;
w_3 = 2 * pi * f_3;
w_4 = 2 * pi * f_4;

a_min = 27 + a_3 * (5/9);
a_max = 0.4 + a_4 / 36;
W_p = 1;
W_s = (w_2-w_1)/(w_4-w_3);
bw = w_2 - w_1;

n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(W_s));
n = ceil(n);

e = sqrt(10^(a_max/10)-1);
a = (1/n)*asinh(1/e);
W_hp = cosh(acosh(1/e)/n);
%% Butterworth angles
y_k1 = 22.5;
y_k2 = -22.5;
y_k3 = 67.5;
y_k4 = -67.5;

%% Chebyshev poles
s_1 = -sinh(a)*cosd(y_k1) + 1i*cosh(a)*sind(y_k1);
s_2 = -sinh(a)*cosd(y_k2) + 1i*cosh(a)*sind(y_k2);
s_3 = -sinh(a)*cosd(y_k3) + 1i*cosh(a)*sind(y_k3);
s_4 = -sinh(a)*cosd(y_k4) + 1i*cosh(a)*sind(y_k4);

W_01 = abs(s_1);
Q1 = W_01/(2*abs(real(s_1)));
W_01 = 1 / W_01;
y_1 = acosd(1/(2*W_01));

W_02 = abs(s_3);
Q2 = W_02/(2*abs(real(s_3)));
W_02 = 1 / W_02;
y_2 = acosd(1/(2*W_02));

%% New poles
p_1 = -W_01/(2*Q1)+j*sqrt(W_01^2 - (W_01/(2*Q1))^2);
p_3 = -W_02/(2*Q2)+j*sqrt(W_02^2 - (W_02/(2*Q2))^2);

%% first pole transformation
S_2 = abs(real(p_1));
W_2 = abs(imag(p_1));
% Geffe algorithm
C = S_2^2 + W_2^2;
q_c = w_0 / bw;
D = 2 * S_2 / q_c;
E = 4 + C / q_c^2;
G = sqrt(E^2 - 4*D^2);
Q_1 = 1 / D * sqrt((E+G)/2);
k = S_2 * Q_1 / q_c;
W = k + sqrt(k^2-1);
w_02 = W * w_0;
w_01 = w_0 / W;

%% second pole transformation
S_2 = abs(real(p_3));
W_2 = abs(imag(p_3));
% Geffe algorithm
C = S_2^2 + W_2^2;
q_c = w_0 / bw;
D = 2 * S_2 / q_c;
E = 4 + C / q_c^2;
G = sqrt(E^2 - 4*D^2);
Q_2 = 1 / D * sqrt((E+G)/2);
k = S_2 * Q_2 / q_c;
W = k + sqrt(k^2-1);
w_04 = W * w_0;
w_03 = w_0 / W;

%%
% we use BoctorHighPass circuits
% we need capacitor 0.01 ìF
% we need gain 0 dB
%% First unit -LPN
w_z1 = w_0 / w_01;
%0.5581<k_11<1
k_11=0.7;
%C1 = 1 / (2 * Q_1);
R11 = 2/(k_11*w_z1^2-1);
R12 = 1 / (1 -k_11);
R13 = (k_11/(                                       Q_1)^2+k_11*w_z1^2-1)/2;
R14 = 1 / k_11;
R15 = 1;
R16 = 1;
C11 = k_11 / (2*Q_1) ;
C12 = 2*Q_1;

%k_12 = 1/(R13+1);
k_12 = 2 / (k_11/Q_1^2+k_11*w_z1^2+1);
k_f1 = w_01;
k_m1 = 10^8 * C11 / k_f1;

C11 = 10^(-8);
R11 = R11 * k_m1;
R12 = R12 * k_m1;
R13 = R13 * k_m1;
R14 = R14 * k_m1;
R15 = R15 * k_m1;
R16 = R16 * k_m1;
C12 = C12/(k_m1 * k_f1);

%k1 = k_12 * (w_0)^2;
k1 = k_12 * (w_z1)^2;

%% Second Unit-HPN
% you should use a HighPassNotch
w_z2 = w_0 / w_02;
k_21 = 1 / w_z2^2 - 1 ;
k_22 = ((2 + k_21)*Q_1^2)/(Q_1^2*(k_21 + 2) + 1);
C2 = 1 / (Q_1*(2 + k_21));
C23 = k_21 * C2;
R21 = 1;
R22 = Q_1^2 * (k_21 + 2)^2;
R23 = 1;
R24 = Q_1^2 * (k_21 + 2);

k_f2 = w_02;
k_m2 = 10^8 * C2 / k_f2;

C2 = 10^(-8);
R21 = R21 * k_m2;
R22 = R22 * k_m2;
R23 = R23 * k_m2;
R24 = R24 * k_m2;
C23 = C23 /(k_m2 * k_f2);

%k2 = k_22 * (w_0)^2;
k2 = k_22 * (w_z2)^2;

%% Third unit -LPN
w_z3 = w_0 / w_03;
%0.7429<k_31<1
k_31=0.8;
C31 = k_31 / (2 * Q_2);
R31 = 2/(k_31*w_z3^2-1);
R32 = 1 / (1 -k_31);
R33 = (k_31/(Q_2)^2+k_31*w_z3^2-1)/2;
R34 = 1 / k_31;
R35 = 1;
R36 = 1;
C31 = k_31 / (2*Q_2) ;
C32 = 2*Q_2;

k_32 = 2/(k_31/Q_2^2 + k_31*w_z3^2+1);

k_f3 = w_03;
k_m3 = 10^8 * C31 / k_f3;

C31 = 10^(-8);
R31 = R31 * k_m3;
R32 = R32 * k_m3;
R33 = R33 * k_m3;
R34 = R34 * k_m3;
R35 = R35 * k_m3;
R36 = R36 * k_m3;
C32 = C32/(k_m3 * k_f3);

%k3 = k_32* (w_0)^2;
k3 = k_32* (w_z3)^2;
%% Forth unit -HPN
% you should use a HighPassNotch
w_z4 = w_0 / w_04;
k_41 = 1/w_z4^2 - 1 ;
k_42 = ((2 + k_41)*Q_2^2)/(Q_2^2*(k_41 + 2) + 1);
C4 = 1 /(Q_2*(2 + k_41));
R41 = 1;
R42 = Q_2^2 * (k_41 + 2)^2;
R43 = 1;
R44 = Q_2^2 * (k_41 + 2);
C42 = k_41 * C4; 

k_f4 = w_04;
k_m4 = 10^8 * C4 / k_f4;

C4 = 10^(-8);
R41 = R41 * k_m4;
R42 = R42 * k_m4;
R43 = R43 * k_m4;
R44 = R44 * k_m4;
C42 = C42 /(k_m4 * k_f4);

%k4 = k_42* (w_0)^2;
k4 = k_42 * (w_z4)^2;

%% Transfer Function
T1 = tf([k_12 0 k_12*w_0^2],[1 w_01/Q_1 w_01^2]);

T2 = tf([k2 0 k2*w_0^2],[1 w_02/Q_1 w_02^2]);
 
T3 = tf([k_32 0 k_32*w_0^2],[1 w_03/Q_2 w_03^2]);

T4 = tf([k4 0 k4*w_0^2],[1 w_04/Q_2 w_04^2]);

T_all = T1*T2*T3*T4;

%% Total Gain
K_tot = k_12*k2*k_32*k4;
gain = abs(evalfr(T_all, w_0 * 1i));
a_tf = (1)/K_tot;
% a>1
R_a = 10000/(a_tf-1);
R_b = 10000;

T_total = a_tf*T_all;

plot_transfer_function(T1, [f_1 f_2 f_3 f_4])

plot_transfer_function(T2, [f_1 f_2 f_3 f_4])

plot_transfer_function(T3, [f_1 f_2 f_3 f_4])

plot_transfer_function(T4, [f_1 f_2 f_3 f_4])

plot_transfer_function(T_total, [f_1 f_2 f_3 f_4])
 
ltiview({'bodemag'}, T1, T2, T3, T4,T_total)

InvSys = inv (T_all);
plot_transfer_function(InvSys, [f_1 f_2 f_3 f_4])
 
InvSys_new1 = inv (T_total);
plot_transfer_function(InvSys_new1, [f_0 f_1 f_2 f_3 f_4])

%% Fourier Analysis
% a_4=9
f11 = (w_0- (w_0-w_3)/2)/(2*pi);
f12 = (w_0+(w_0+w_3)/3)/(2*pi);
f13 = 0.4*w_1 /(2*pi);
f14 = (2.5*w_2) / (2*pi);
f15 = (3*w_2)/(2*pi) ;
fs = 200*10^3;
T = 0.002;
dt = 1/fs;
t = 0:dt:(T);

u1 = 0.5*cos(2*pi*f11*t)+0.8* cos(2*pi* f12*t)+0.8*cos(2*pi* f13*t)+0.6*cos(2*pi*f14*t)+1.2*cos(2*pi*f15*t);
figure
plot(u1)

N=T/dt;
figure
lsim(T_total,u1,t)
xt=lsim(T_total,u1,t);
figure
plot(t,xt)
n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
figure
plot(f,p1)

nfft=n;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*fs/nfft; 
figure
plot(f,y_mag)



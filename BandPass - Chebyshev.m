%% Tzikas Tryfon Rigas
%% 8589

%% BandPass Chebyshev


clear
clc

%initializations
a_1 = 8;
a_2 = 5;
a_3 = 8;
a_4 = 9;

f_0 = 900;
f_1 = 650 + 25 * a_3;
f_2 = f_0^2 / f_1;
D = 2.2 * (f_0^2-f_1^2) / f_1;
f_3 = (-D+sqrt(D^2+4*f_0^2))/2;
f_4 = f_0^2/f_3;

w_0 = 2 * pi * f_0;
w_1 = 2 * pi * f_1;
w_2 = 2 * pi * f_2;
w_3 = 2 * pi * f_3;
w_4 = 2 * pi * f_4;

W_p = 1;
W_s = (w_4-w_3)/(w_2-w_1);
bw = w_2 - w_1;
a_min = 28.5 + a_4 * (5/9);
a_max = 0.5 + a_3 / 36;

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
s_1 = -sinh(a)*cosd(y_k1) + j*cosh(a)*sind(y_k1);
s_2 = -sinh(a)*cosd(y_k2) + j*cosh(a)*sind(y_k2);
s_3 = -sinh(a)*cosd(y_k3) + j*cosh(a)*sind(y_k3);
s_4 = -sinh(a)*cosd(y_k4) + j*cosh(a)*sind(y_k4);


%% Metasxhmatismos prwtou polou
S_2 = abs(real(s_1));
W_2 = abs(imag(s_1));
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

%% Metasxhmatismos deuterou polou
S_2 = abs(real(s_3));
W_2 = abs(imag(s_3));
% Geffe algorithm
C = S_2^2 + W_2^2;
q_c = w_0 / bw;
D = 2 * S_2 / q_c;
E = 4 + C / q_c^2;
G = sqrt(E^2 - 4*D^2);
Q_2 = sqrt((E+G)/2)/D;
k = S_2 * Q_2 / q_c;
W = k + sqrt(k^2-1);
w_04 = W * w_0;
w_03 = w_0 / W;

%%
% a3=8 ara xrhsimopoioume strathgikh 1
% a2=5 ara theloume puknwth 1.0 ìF
% a4=9 ara kerdos 5 dB
% omws Q>>5 ara tha xrhsimopoihsoyme Q enchansement
%% Prwth Monada
C11 = 1;
C12 = 1;
b = 50;
R11 = 1 / sqrt(b);
R12 = sqrt(b);
k = (Q_1*(b+2)-sqrt(b))/(2*Q_1-sqrt(b));
R1A = 1000;
R1B = (k - 1) * R1A;
H1 = k*b/(2*(b-1)-b);

k_f = w_01;
k_m = 10^6 * C11 / k_f; 

C11 = 10^(-6);
C12 = 10^(-6);
R11 = R11 * k_m;
R12 = R12 * k_m;
R1A = R1A * k_m;
R1B = R1B * k_m;
% we want gain: 5dB so 20*log(X) = 5 
%k_desired = 10^0.25;
% tha prepei na ginei exasthenisi tou kerdous
%a = k_desired / H1;

%Z11 = R11 / a;
%Z12 = R12 /(1-a);

%% Deuterh Monada
C21 = 1;
C22 = 1;
b = 50;
R21 = 1 / sqrt(b);
R22 = sqrt(b);
k = (Q_1*(b+2)-sqrt(b))/(2*Q_1-sqrt(b));
R2A = 1000;
R2B = (k - 1) * R2A;
H2 = k*b/(2*(b-1)-b);

k_f = w_02;
k_m = 10^6 * C21 / k_f; 

C21 = 10^(-6);
C22 = 10^(-6);
R21 = R21 * k_m;
R22 = R22 * k_m;
R2A = R2A * k_m;
R2B = R2B * k_m;
% we want gain: 5dB so 20*log(X) = 5 
%k_desired = 1;
% tha prepei na ginei exasthenisi tou kerdous
%a = k_desired / H2;

%Z21 = R21 / a;
%Z22 = R22 /(1-a);

%% Trith Monada
C31 = 1;
C32 = 1;
b = 50;
R31 = 1 / sqrt(b);
R32 = sqrt(b);
k = (Q_2*(b+2)-sqrt(b))/(2*Q_2-sqrt(b));
R3A = 1000;
R3B = (k - 1) * R3A;
H3 = k*b/(2*(b-1)-b);

k_f = w_03;
k_m = 10^6 * C31 / k_f; 

C31 = 10^(-6);
C32 = 10^(-6);
R31 = R31 * k_m;
R32 = R32 * k_m;
R3A = R3A * k_m;
R3B = R3B * k_m;
% we want gain: 5dB so 20*log(X) = 5 
%k_desired = 1;
% tha prepei na ginei exasthenisi tou kerdous
%a = k_desired / H3;

%Z31 = R31 / a;
%Z32 = R32 /(1-a);

%% Tetarth Monada
C41 = 1;
C42 = 1;
b = 50;
R41 = 1 / sqrt(b);
R42 = sqrt(b);
k = (Q_2*(b+2)-sqrt(b))/(2*Q_2-sqrt(b));
R4A = 1000;
R4B = (k - 1) * R4A;
H4 = k*b/(2*(b-1)-b);

k_f = w_04;
k_m = 10^6 * C41 / k_f; 

C41 = 10^(-6);
C42 = 10^(-6);
R41 = R41 * k_m;
R42 = R42 * k_m;
R4A = R4A * k_m;
R4B = R4B * k_m;

% we want gain: 5dB so 20*log(X) = 5 
%k_desired = 1;
% tha prepei na ginei exasthenisi tou kerdous
%a = k_desired / H4;

%Z41 = R41 / a;
%Z42 = R42 /(1-a);

%% Transfer Function
T1 = tf( [0 H1*(w_01/Q_1) 0], [ 1, ( w_01 / Q_1), w_01^2 ] );

T2 = tf( [0 H2*(w_02/Q_1) 0], [ 1, ( w_02 / Q_1), w_02^2 ] );

T3 = tf( [0 H3*(w_03/Q_2) 0], [ 1, ( w_03 / Q_2), w_03^2 ] );

T4 = tf( [0 H4*(w_04/Q_2) 0], [ 1, ( w_04 / Q_2), w_04^2 ] );

T_all = T1*T2*T3*T4;

%% Total Gain
K_total = H1*H2*H3*H4;
gain = abs(evalfr(T_all, w_0 * 1i));

a = 10^0.25/gain;
a1 = 10^0.25/K_total;
% a<1
Z1 = R11/(10^0.25/K_total);
Z2 = R11/(1-((10^0.25)/K_total));

T_total = a * T_all;

plot_transfer_function(T1, [f_1 f_2 f_3 f_4])

plot_transfer_function(T2, [f_1 f_2 f_3 f_4])

plot_transfer_function(T3, [f_1 f_2 f_3 f_4])

plot_transfer_function(T4, [f_1 f_2 f_3 f_4])

plot_transfer_function(T_total, [f_1 f_2 f_3 f_4])

ltiview({'bodemag'}, T1, T2, T3, T4, T_total)

InvSys_new = inv (T_all);
plot_transfer_function(InvSys_new, [f_1 f_2 f_3 f_4])

InvSys_new1 = inv (T_total);
plot_transfer_function(InvSys_new1, [f_1 f_2 f_3 f_4])

%% Fourier Analysis
% a_4 = 9
f11 = (w_0 - (w_0-w_1)/2) / (2*pi);
f12 = (w_0 + (w_0+w_1)/3) / (2*pi);
f13 = 0.4*w_3 / (2*pi);
f14 = (2.5*w_4) / (2*pi);
f15 = (3*w_4) / (2*pi) ;
fs= 200*10^3;
T=0.002;
dt=1/fs;
t=0:dt:(T);

u1 = cos(2*pi*f11*t)+0.8* cos(2*pi* f12*t)+0.8*cos(2*pi* f13*t)+0.6*cos(2*pi*f14*t)+0.5*cos(2*pi*f15*t);
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
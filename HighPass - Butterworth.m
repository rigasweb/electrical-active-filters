% Tzikas Tryfon Rigas

% HighPass Butterworth


clear

%initializations
a_1 = 8;
a_2 = 5;
a_3 = 8;
a_4 = 9;

m = 2;

f_p = ( 4 + m ) * 10^3 ;
f_s = f_p / 2.6 ;
a_min = 24 + a_3 * (6/9);
a_max = 0.5 + a_4 / 36;


w_p = 2 * pi * f_p;
w_s = 2 * pi * f_s;

W_p = 1;
W_s = w_p / w_s;

e = sqrt(10^(a_max/10)-1);
n = (log10( (10^(a_min/10)-1) / (10^(a_max/10)-1) ) / (2*log10(W_s/W_p)) );

n = ceil(n);

W_0= W_p / (10^(a_max/10) - 1)^(1/(2*n));
w_0 = w_p / W_0;
f_0 = w_0 / (2*pi);

% butterworth angles
y_k1 = 0;
y_k2 = +36;
y_k3 = -36;
y_k4 = +72;
y_k5 = -72;

% poloi
p_1 = -1;
p_2 = -0.809 + j*0.5877;
p_3 = -0.809 - j*0.5877;
p_4 = -0.309 + j*0.9510;
p_5 = -0.309 - j*0.9510;

% w0, Q of poles 
Q_1 = 0.5;

Q_23 = 0.618;

Q_45 = 1.618;

% we need capacitor 0.01 ìF
% we need gain 10 dB
% First unit 
C11 = 1; 
R11 = 1;
k_f1 = w_0;
p1 = 1;
C11 = 0.01 * 10^(-6);
k_m1 = 1 / (k_f1 * C11);
R11 = R11 * k_m1;

% Second Unit
C21 = 1;
C22 = 1;
R21 = 1/(2*Q_23);
R22 = 2*Q_23;
C22 = 0.01 * 10^(-6);
C21 = C22;
k2 = 1;
k_f2 = w_0;
k_m2 = 1 / (k_f2 * C22);
R21 = R21 * k_m2;
R22 = R22 * k_m2;

% Third unit
C31 = 1;
C32 = 1;
R31 = 1 / (2*Q_45);
R32 = 2*Q_45;
k3 = 1;
k_f3 = w_0;
C32 = 0.01 * 10^(-6);
C31 = C32;
k_m3 = 1 / (k_f3 * C32);
R31 = k_m3 * R31;
R32 = k_m3 * R32;

% Transfer functions
T1 = tf ([1 0], [1 w_0]);

T2 = tf ([ k2 0 0], [1 w_0/Q_23 w_0^2]);

T3 = tf ([ k3 0 0], [1 w_0/Q_45 w_0^2]);

T_all = T1*T2*T3 ;

% Total Gain
gain = abs(evalfr(T_all, w_0 * 1i));
a_tf = (10^0.5)/gain;
% a>1
R_b = 1 * 10^3;
R_a = R_b*(a_tf-1);


T_total = a_tf * T_all;

plot_transfer_function(T1, [f_p, f_s])

plot_transfer_function(T2, [f_p, f_s])

plot_transfer_function(T3, [f_p, f_s])

plot_transfer_function(T_total, [f_p, f_s])

ltiview({'bodemag'}, T1, T2, T3, T_total)
 
InvSys_new = inv (T_all)
plot_transfer_function(InvSys_new, [f_p f_s])

InvSys_new1 = inv (T_total)
plot_transfer_function(InvSys_new1, [f_p f_s])

% Fourier Analysis
f11= (0.2*w_s) / (2*pi);
f12= (0.7*w_s) / (2*pi);
f13= (1.6*w_p) / (2*pi);
f14= (2.4*w_p) / (2*pi);
f15= (3.5*w_p) / (2*pi) ;
fs= 200*10^3;
T=0.002;
dt=1/fs;
t=0:dt:(T);

u1= cos(2*pi*f11*t)+0.6*cos(2*pi*f12*t)+1.5*cos(2*pi*f13*t)+0.7*cos(2*pi*f14*t)+0.4*cos(2*pi*f15*t);
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




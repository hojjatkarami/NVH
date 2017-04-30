function dx = eng_mount(t, x, M, C, K)

global T_Amp T_Freq T_phase

dx = x;
dx(1:6) = x(7:12);

T = [0; 0; 0; 0;T_Amp*sin(T_Freq*t+T_phase);0];      % Exerted torque by the crankshaft

dx(7:12) = M\(-K*x(1:6)-C*x(7:12)+T);
% for now this is not a class... just a bunch of defintions !

clear vehicle;

vehicle.dt = 1/100;
vehicle.m = 4.05;
vehicle.weight = vehicle.m*9.81;
vehicle.tmin = 0;
vehicle.tmax = vehicle.weight * 1.2 / 4;
vehicle.Iyy = .542;
vehicle.motors_fusloc = [0.28575,   0.50165   , 0.71755   ,    0.93345];
vehicle.cg = .6096;
vehicle.motors_xloc = -vehicle.motors_fusloc + vehicle.cg;
%[Fx Fz My] = LCM * t
vehicle.lcm = [0 0 0 0 ; -1 -1 -1 -1 ; vehicle.motors_xloc / vehicle.Iyy]; 
vehicle.FM = [0 ; 0 ; 0];
vehicle.U = zeros(4,1);

vehicle.x = [0; 0]; %theta , theta_dot
vehicle.A = [0 1; 0 0];
vehicle.B = [zeros(1,4) ; vehicle.lcm(3,:)];
vehicle.sysd = c2d(ss(vehicle.A,vehicle.B,[],[]),vehicle.dt);

 %command and command rate for filter
vehicle.control_pid.theta_cmd = 0;
vehicle.control_pid.theta_cmd_rate = 0;
vehicle.control_pid.theta_errI = 0; %integrated error 
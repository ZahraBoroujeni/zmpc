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
vehicle.lcm = [0 0 0 0 ; -1 -1 -1 -1 ; vehicle.motors_xloc]; 
vehicle.FM = [0 ; 0 ; 0];
vehicle.U = zeros(4,1);


% Gy_reduced = plantReduce(Gy,{'\theta','q','u'},'My_{cmd}',{'\theta','q_{INS}','u'});
% vehicle.A =  Gy_reduced.a;
vehicle.A = [0 1; 0 0];
vehicle.B = [0 ; 1 / vehicle.Iyy];
vehicle.sysd = c2d(ss(vehicle.A,vehicle.B,eye(size(vehicle.A)),[]),vehicle.dt);
vehicle.sysd.InputName = 'My';
vehicle.sysd.OutputName = {'\theta','q'};
 
Glcm = ss(vehicle.lcm,'InputName',{'t1','t2','t3','t4'}, ...
                        'OutputName',{'Fx_{cmd}','Fz_{cmd}','My_{cmd}'});

                    
          
%Motor dynamics:
M_omega1 = 22/1.4;
M_omega2 = 60*2*pi;
motor_dynamics = tf(1/( (s/M_omega1) + 1) /( (s/M_omega2) + 1));
motor_dynamics = tf(1);
motor_dynamics = ss(motor_dynamics);
z = tf('z',.01);
%motor_dynamics = c2d(motor_dynamics,.01) * ss(1/(z^2)); %.02 time delay;
motor_dynamics.InputName = 'My_{cmd}';
motor_dynamics.OutputName = 'My';

vehicle.actD.sysd = connect(motor_dynamics, Glcm, Glcm.InputName, motor_dynamics.OutputName);
vehicle.actD.tstates = zeros(size(vehicle.actD.sysd.a,1),1);

vehicle.sysd = connect(vehicle.actD.sysd,vehicle.sysd,vehicle.actD.sysd.InputName, vehicle.sysd.OutputName);

%converge the actuator to trim state initially!
tin = vehicle.weight/4 * ones(4,1);
for i = 1:100
    vehicle.actD.tstates = vehicle.actD.sysd.a * vehicle.actD.tstates + vehicle.actD.sysd.b * tin;
end
          
                    
             
                    
vehicle.sysSim = connect(Gy,Gz,Glcm, [Glcm.InputName; 'My_{dist}'; 'u_{wind}'],[Gy.OutputName; Gz.OutputName]); 
vehicle.sysSim = c2d(vehicle.sysSim,vehicle.dt);
vehicle.xsim = 0 * zeros(size(vehicle.sysSim.a,1),1);

vehicle.thetaIndex = find(strcmp(vehicle.sysSim.OutputName,'\theta'));
vehicle.qIndex = find(strcmp(vehicle.sysSim.OutputName,'q_{INS}'));
vehicle.uIndex = find(strcmp(vehicle.sysSim.OutputName,'u'));
vehicle.Myindex = find(strcmp(vehicle.sysSim.OutputName,'My'));
vehicle.Fxindex = find(strcmp(vehicle.sysSim.OutputName,'Fx}'));
vehicle.Fzindex = find(strcmp(vehicle.sysSim.OutputName,'Fz'));


vehicle.Fzdist = 0 ;
vehicle.Mydist = 0 ;
vehicle.uwind = 0;



%same as above but bunching actuators as My
vehicle.AMy = vehicle.A(1:2,1:2);
vehicle.BMy = [ 0 ; 1/vehicle.Iyy];
vehicle.sysdMy = c2d(ss(vehicle.AMy,vehicle.BMy,[],[]),vehicle.dt);

%command and command rate for filter
vehicle.control_pid.theta_cmd_wn = 18;
vehicle.control_pid.theta_cmd_z = 1;
vehicle.control_pid.theta_cmd = 0;
vehicle.control_pid.theta_cmd_rate = 0;
vehicle.control_pid.theta_errI = 0; %integrated error 

vehicle.control_lqru = vehicle.control_pid;
%this is to handle the fact that we are doing tracking instead of
%regulating
vehicle.control_lqru.BinvIA = pinv(vehicle.sysdMy.B)*(eye(size(vehicle.sysdMy.A)) - vehicle.sysdMy.A);
vehicle.control_lqru.Q = diag([1 .01]);
vehicle.control_lqru.R = .01;
%vehicle.control_lqru.K = dlqr(vehicle.sysdMy.A,vehicle.sysdMy.B,...
 %                            vehicle.control_lqru.Q, vehicle.control_lqru.R);

%Testing finite horizon lqr...
vehicle.control_lqru.K = lqr_finite(vehicle.sysdMy.A,vehicle.sysdMy.B,300,...
                             vehicle.control_lqru.Q, vehicle.control_lqru.R);

%control mix weight on Fz_err and My_err
vehicle.ctrlmix_W = [1 0; 0 1000];


vehicle.control_cvx.solved = false;
vehicle.control_cvx.Ftarget = [0; 0];


%Estimator settings
vehicle.estimator_dist.Myd = 0; 
vehicle.estimator_dist.xd = [0 ; 0 ;vehicle.estimator_dist.Myd]; %theta, theta_dot, Myd

vehicle.estimator_dist.A = [vehicle.sysdMy.A , vehicle.sysdMy.B; zeros(1,2) , 1];
vehicle.estimator_dist.B = [vehicle.sysdMy.B ; 0];

vehicle.estimator_dist.Cxd = [eye(2) , zeros(2,1)];
vehicle.estimator_dist.Lxd = -[1 0; 
                              0 1;
                              0 5];


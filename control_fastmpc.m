function vehicle = control_fastmpc(vehicle, Ftarget_in)

%x = Ax * Bu where x = [theta, theta_dot]
%y = [Fx; delta_Fz] = [-W*theta ; ]
%y = [-W , 0] * x + [-1 -1 -1 -1] * u
T = 100;
n = size(vehicle.sysd.a,1);
m = size(vehicle.sysd.b,2);

%initialize the state vector if it's not already initialized
if(~isfield(vehicle.control_cvx,'x0'));
    vehicle.control_cvx.x0  = zeros(n,1);
    vehicle.control_cvx.x0(2) = -vehicle.weight/vehicle.sysd.c(2,2);
    vehicle.control_cvx.thetacmd = 0;
    vehicle.control_cvx.cmd_rate = 0;
end
%overwrite it with the current measured values
x0  = vehicle.control_cvx.x0;
x0(vehicle.sysdthetaIndex) = vehicle.theta;
x0(vehicle.sysdqIndex) = vehicle.q;
x0(vehicle.sysduIndex) = vehicle.u;

%define the costs
Qy = diag([5 ; 50]);
Qyf = diag([5; 50]);
r =  1;
R = diag(r*ones(m,1));
C = vehicle.sysd.C;
C(2,:) = C(2,:) * cos(vehicle.theta);
D = vehicle.sysd.D;

    
if(~isfield(vehicle.control_cvx,'usol'))
    vehicle.control_cvx.usol = ones(m,T) * vehicle.weight/4;
    vehicle.control_cvx.xsol = repmat(x0,1,T);
    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.numSol = 0;
    vehicle.control_cvx.deltaT = 0;
    vehicle.control_cvx.dtHist = [];
    vehicle.control_cvx.solved = false;
else
   
   delta = vehicle.control_cvx.iter + 1;
   %shift the solution
   vehicle.control_cvx.usol(:,1:(end-delta+1)) = vehicle.control_cvx.usol(:,delta:end);
   vehicle.control_cvx.xsol(:,1:(end-delta+1)) = vehicle.control_cvx.xsol(:,delta:end);
end


if( 1 || (~vehicle.control_cvx.solved))

    
    %compute the noise
    w = zeros(n,1);
    wdelta = vehicle.sysdMy.b *  vehicle.estimator_dist.Myd;
    w(vehicle.sysdthetaIndex) = wdelta(1);
    w(vehicle.sysdqIndex) = wdelta(2);
    
    
    %compute the desired y
    Ftarget = Ftarget_in;
    wn = 18;
    zeta = 1;
    thetacmd = -asin(Ftarget(1)/vehicle.weight);
    thetamax = 15*pi/180;
    thetacmd = max(min(thetacmd , thetamax),-thetamax);
    cmd_accel = wn^2 * (thetacmd - vehicle.control_cvx.thetacmd)...
        - 2.0 * zeta * wn *vehicle.control_cvx.cmd_rate ;
    vehicle.control_cvx.cmd_rate = vehicle.control_cvx.cmd_rate + ...
        cmd_accel * vehicle.dt;
    vehicle.control_cvx.cmd_rate = clip(vehicle.control_cvx.cmd_rate,-vehicle.control_pid.slewrate, vehicle.control_pid.slewrate);
    
    vehicle.control_cvx.thetacmd = vehicle.control_cvx.thetacmd + ...
        vehicle.control_cvx.cmd_rate * vehicle.dt;
    
    Ftarget(1) = -sin(vehicle.control_cvx.thetacmd)*vehicle.weight;

    ydesired = Ftarget;
    
    
    
%    function [U, X] = fastmpcsolve(A,B,C,D,Qy,Qyf,R, yd, w, x0, X0, U0, ubounds, xbounds)
    tstart = tic;

    [U, X] = fastmpcsolve(vehicle.sysd.a, vehicle.sysd.b, C, D, ...
                          Qy, Qyf, R, T,...
                          ydesired, w, x0,...
                          vehicle.control_cvx.usol, vehicle.control_cvx.xsol, ...
                          [vehicle.tmin, vehicle.tmax], []);
                      
    
    dt = toc(tstart);
    vehicle.control_cvx.dtHist = [    vehicle.control_cvx.dtHist ; dt];
    vehicle.control_cvx.deltaT = vehicle.control_cvx.deltaT +dt ;
    vehicle.control_cvx.numSol = vehicle.control_cvx.numSol + 1;
    
    vehicle.control_cvx.usol = U;
    vehicle.control_cvx.xsol = X;
    vehicle.control_cvx.F = C* vehicle.control_cvx.xsol + D * vehicle.control_cvx.usol;
    vehicle.control_cvx.Ftarget = Ftarget_in;

    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.solved = true;
end


vehicle.control_cvx.iter = vehicle.control_cvx.iter + 1;
if(vehicle.control_cvx.iter < size(vehicle.control_cvx.usol,2))
    vehicle.U = vehicle.control_cvx.usol(:,vehicle.control_cvx.iter);
    vehicle.control_cvx.x0 = vehicle.sysd.a *  x0  +  vehicle.sysd.b * vehicle.U;
end

end
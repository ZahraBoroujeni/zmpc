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
C(2,:) = C(2,:) / cos(vehicle.theta);
D = vehicle.sysd.D;

    
if(~isfield(vehicle.control_cvx,'usol'))
    vehicle.control_cvx.usol = ones(m,T) * vehicle.weight/4;
    vehicle.control_cvx.xsol = repmat(x0,1,T);
    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.numSol = 0;
    vehicle.control_cvx.deltaT = 0;
    vehicle.control_cvx.solved = false;
else
   
   delta = vehicle.control_cvx.iter + 1;
   %shift the solution
   vehicle.control_cvx.usol(:,1:(end-delta+1)) = vehicle.usol(:,delta:end);
   vehicle.control_cvx.xsol(:,1:(end-delta+1)) = vehicle.xsol(:,delta:end);
end


if( 1 || vehicle.control_cvx.iter == 2 || (~vehicle.control_cvx.solved))

    
    %compute the noise
    w = zeros(n,1);
    wdelta = vehicle.sysdMy.b *  vehicle.estimator_dist.Myd;
    w(vehicle.sysdthetaIndex) = wdelta(1);
    w(vehicle.sysdqIndex) = wdelta(2);
    
    
    %compute the desired y
    Ftarget = Ftarget_in;
    theta_Fx = asind(Ftarget(1)/vehicle.weight);
    theta_Fx = max(min(theta_Fx,15),-15);
    Ftarget(1) = sind(theta_Fx)*vehicle.weight;    
    ydesired = Ftarget;
    
    
%    function [U, X] = fastmpcsolve(A,B,C,D,Qy,Qyf,R, yd, w, x0, X0, U0, ubounds, xbounds)
    tstart = tic;

    [U, X] = fastmpcsolve(sysd.a, sysd.b, C, D, ...
                          Qy, Qyf, R, ydesired, w, x0,...
                          vehicle.control_cvx.xsol, vehicle.control_cvx.usol, ...
                          [vehicle.tmin, vehicle.tmax], []);
                      
    
    vehicle.control_cvx.deltaT = vehicle.control_cvx.deltaT + toc(tstart);
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
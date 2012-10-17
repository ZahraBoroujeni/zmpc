function vehicle = control_lqru(vehicle, Ftarget)
wn = vehicle.control_lqru.theta_cmd_wn;
zeta = vehicle.control_lqru.theta_cmd_z;

K = vehicle.control_lqru.K;
BinvIA = vehicle.control_lqru.BinvIA;

thetamax = 15*pi/180;
thetamin = -thetamax;

theta = vehicle.x(1);
q = vehicle.x(2);

Fxcmd = Ftarget(1);
Fzcmd = Ftarget(2)/cos(theta); %Fz / cos(theta)

thetacmd = -Fxcmd/vehicle.weight;
thetacmd = max(min(thetacmd , thetamax),thetamin);


cmd_accel = wn^2 * (thetacmd - vehicle.control_lqru.theta_cmd)...
             - 2.0 * zeta * wn *vehicle.control_lqru.theta_cmd_rate ;
         
vehicle.control_lqru.theta_cmd_rate = vehicle.control_lqru.theta_cmd_rate + ...
              cmd_accel * vehicle.dt;
          
vehicle.control_lqru.theta_cmd = vehicle.control_lqru.theta_cmd + ...
    vehicle.control_lqru.theta_cmd_rate * vehicle.dt;

xcmd = [vehicle.control_lqru.theta_cmd ; vehicle.control_lqru.theta_cmd_rate];
xerr = (xcmd - vehicle.x);

vehicle.control_lqru.Mycmd = K * xerr - BinvIA * xcmd;
vehicle.U = control_mix(vehicle, [Fzcmd ;  vehicle.control_lqru.Mycmd]);
%vehicle.U = min(max(vehicle.U,vehicle.tmin),vehicle.tmax);

%backsolving ...


end
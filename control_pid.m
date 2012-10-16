function vehicle = control_pid(vehicle, Ftarget)
kp = 9.55;
kd = 2.743;
kd2 = 1;
ki = 4.775;
wn = 18; zeta = 1;

thetamax = 15*pi/180;
thetamin = -thetamax;

theta = vehicle.x(1);
q = vehicle.x(2);

Fxcmd = Ftarget(1);
Fzcmd = Ftarget(2)/cos(theta); %Fz / cos(theta)

thetacmd = -Fxcmd/vehicle.weight;
thetacmd = max(min(thetacmd , thetamax),thetamin);


cmd_accel = wn^2 * (thetacmd - vehicle.control_pid.theta_cmd)...
             - 2.0 * zeta * wn *vehicle.control_pid.theta_cmd_rate ;
         
vehicle.control_pid.theta_cmd_rate = vehicle.control_pid.theta_cmd_rate + ...
              cmd_accel * vehicle.dt;
          
vehicle.control_pid.theta_cmd = vehicle.control_pid.theta_cmd + ...
    vehicle.control_pid.theta_cmd_rate * vehicle.dt;


err = (vehicle.control_pid.theta_cmd - theta);
vehicle.control_pid.theta_errI = vehicle.control_pid.theta_errI + err*vehicle.dt;

vehicle.control_pid.Mycmd = kp * err + ...
        kd * (vehicle.control_pid.theta_cmd_rate - q) ...
        - kd2 * q ...
        + ki * vehicle.control_pid.theta_errI;

vehicle.U = pinv(vehicle.lcm) * [0 ; Fzcmd ;  vehicle.control_pid.Mycmd];
%vehicle.U = min(max(vehicle.U,vehicle.tmin),vehicle.tmax);

%backsolving ...


end
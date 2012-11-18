function vehicle = control_u_pid(vehicle, utarget)           

% .OuterVuGains =
%             {
%              /*g q ~ -udotdot*/       0.00f,
%              /*s*u_cmd - udot */      0.00f,
%              /*u_cmd-u*/             -0.16f,
%              /*uErrI*/                0.00f,
%              /*uErrD*/                0.00f,
%             },

kp = -0.16;
wn = vehicle.control_u_pid.cmd_wn;
zeta = vehicle.control_u_pid.cmd_z;

umax = 3;
umin = -3;



ucmd = clip(utarget,umin,umax);
u = vehicle.u;





cmd_accel = wn^2 * (ucmd - vehicle.control_u_pid.cmd)...
             - 2.0 * zeta * wn *vehicle.control_u_pid.cmd_rate ;
         
vehicle.control_u_pid.cmd_rate = vehicle.control_u_pid.cmd_rate + ...
              cmd_accel * vehicle.dt;
          
vehicle.control_u_pid.cmd_rate = clip(vehicle.control_u_pid.cmd_rate,-vehicle.control_u_pid.slewrate, vehicle.control_u_pid.slewrate);           
          
vehicle.control_u_pid.cmd = vehicle.control_u_pid.cmd + ...
    vehicle.control_u_pid.cmd_rate * vehicle.dt;


err = (vehicle.control_u_pid.cmd - u);

vehicle.control_u_pid.FxTarget = - kp * err * vehicle.weight;

end
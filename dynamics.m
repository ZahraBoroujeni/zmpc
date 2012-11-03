function vehicle = dynamics(vehicle,t)

t =  min(max(t,vehicle.tmin),vehicle.tmax);


vehicle.x = vehicle.sysd.A * vehicle.x + vehicle.sysd.B * t + vehicle.sysdMy.B * vehicle.Mydist;

%lcm in body axes then rotated to inertial!
ct = cos(vehicle.x(1));
st = sin(vehicle.x(1));
vehicle.FM = [ct st 0; -st ct 0; 0 0 1] * vehicle.lcm * t; 

end
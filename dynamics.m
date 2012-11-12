function vehicle = dynamics(vehicle,t)

tin =  min(max(t,vehicle.tmin),vehicle.tmax);
%remove trim...
tin = tin - vehicle.weight/4;

uin = [tin ; vehicle.Mydist ; vehicle.uwind];

vehicle.xsim = vehicle.sysSim.A * vehicle.xsim + vehicle.sysSim.B * uin;
vehicle.ysim = vehicle.sysSim.C * vehicle.xsim + vehicle.sysSim.D * uin;


vehicle.theta = vehicle.ysim(vehicle.thetaIndex);
vehicle.q = vehicle.ysim(vehicle.qIndex);
vehicle.u = vehicle.ysim(vehicle.uIndex);
vehicle.My = vehicle.ysim(vehicle.Myindex);
vehicle.bodyFM = [ 0 ; %vehicle.ysim(vehicle.Fxindex);
                   vehicle.ysim(vehicle.Fzindex) - vehicle.weight;
                   vehicle.ysim(vehicle.Myindex)];

%rotate body forces to inertial
ct = cos(vehicle.theta);
st = sin(vehicle.theta);
vehicle.FM = [ct st 0; -st ct 0; 0 0 1] * vehicle.bodyFM; 

end
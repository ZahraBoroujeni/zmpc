function vehicle = estimator_dist(vehicle)

%We assume all of the forces have been simulated and we have measured the
%output

%get the measurement
ym = vehicle.x;
%what's the current expected measurement?
yhat = vehicle.estimator_dist.Cxd * vehicle.estimator_dist.xd;

%My is FM(3)
My = vehicle.FM(3);

vehicle.estimator_dist.xd = vehicle.estimator_dist.A * vehicle.estimator_dist.xd + ...
                           vehicle.estimator_dist.B * My + ...
                           vehicle.estimator_dist.Lxd * (yhat - ym);
                       
                       
vehicle.estimator_dist.Myd = vehicle.estimator_dist.xd(end);                       


end
function vehicle = estimator(vehicle)

%We assume all of the forces have been simulated and we have measured the
%output



ym = [vehicle.theta; vehicle.q];


%what's the current expected measurement?
yhat = vehicle.estimator_dist.Cxd * vehicle.estimator_dist.xd;


vehicle.estimator_dist.xd = vehicle.estimator_dist.A * vehicle.estimator_dist.xd + ...
                           vehicle.estimator_dist.B * vehicle.My + ...
                           vehicle.estimator_dist.Lxd * (yhat - ym);
                       
                       
vehicle.estimator_dist.Myd = vehicle.estimator_dist.xd(end);                       


end
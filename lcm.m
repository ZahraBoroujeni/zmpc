function FM = lcm(vehicle,t)

%clip thrusts to [tmin tmax] then use LCM !

FM = vehicle.lcm * min(max(t,vehicle.tmin),vehicle.tmax);


end
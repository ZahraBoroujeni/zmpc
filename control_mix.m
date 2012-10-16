function U = control_mix(vehicle, FMtarget)

n = 4;
A = vehicle.lcm(2:end,:);
cvx_begin
    variable thrust(n)
    minimize( norm(A*thrust - FMtarget) + .01*norm(thrust)) 
cvx_end


U = thrust;

end
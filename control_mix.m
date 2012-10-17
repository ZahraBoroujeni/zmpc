function U = control_mix(vehicle, FMtarget)
%FMtarget is [Fzcmd ; MyCmd]

n = 4;
A = vehicle.lcm(2:end,:);
W = vehicle.ctrlmix_W;

U = pinv(A) * FMtarget;

%only optimize if we need to...
if(~(all(U >= 0) && all(U <= vehicle.tmax)))
cvx_begin
    variable thrust(n)
    minimize( norm(W*(A*thrust - FMtarget)) + .01*norm(thrust)) 
    subject to
        thrust >= 0
        thrust <= vehicle.tmax
cvx_end
U = thrust;

end

end
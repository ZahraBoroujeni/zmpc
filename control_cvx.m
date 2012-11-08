function vehicle = control_cvx(vehicle, Ftarget_in)


if(0 || (~vehicle.control_cvx.solved) || any(Ftarget_in ~= vehicle.control_cvx.Ftarget))

    Ftarget = Ftarget_in;
    theta_Fx = asind(Ftarget(1)/vehicle.weight);
    theta_Fx = max(min(theta_Fx,15),-15);
    Ftarget(1) = sind(theta_Fx)*vehicle.weight;
    %x = Ax * Bu where x = [theta, theta_dot]
    %y = [Fx; delta_Fz] = [-W*theta ; ]
    %y = [-W , 0] * x + [-1 -1 -1 -1] * u
    N = 100;
    n = 4;
    
    Q = diag([1 10]);
    r = 0.1;
    R = diag(r*ones(n,1));
    
    C = [-vehicle.weight , 0 ; 0, 0];
    D = [0 0 0 0 ; -1 -1 -1 -1];
    
    A = [C'*Q*C , -C'*Q*D ; -D'*Q*C , D'*Q*D + R];
    yd = Ftarget;
    c = -2*[C'*Q; D'*Q] * yd;
    
    Abar = repmat({A},1,N);
    Abar = blkdiag(Abar{:});
    
    cbart = repmat(c,N,1)';
    
    
    cvx_begin
    variable xu(2*N + n*N);
    
    minimize transpose(xu)*Abar*xu + cbart*xu;
    subject to 
        xu(1) == vehicle.x(1);
        xu(2) == vehicle.x(2);
%         xu((N-1)*(2+n)+1) >= -15*pi/180;
%         xu((N-1)*(2+n)+1) <= 15*pi/180;
%         
%         xu((N-1)*(2+n)+2) >= -1*pi/180;
%         xu((N-1)*(2+n)+2) <= 1*pi/180;
        
        for i = 1:N-1
            xi = (i-1)*(2+n)+1;
            xidx = xi:(xi+1);
            xidxp1 = xidx + (2+n);
            ui = xi+2; 
            uidx = ui:(ui+n-1);
            %X(k+1) = A*x + B*u
            xu(xidxp1) == vehicle.sysd.a * xu(xidx) + vehicle.sysd.b * xu(uidx) + vehicle.sysdMy.b * vehicle.estimator_dist.Myd;
            %max thrust
            xu(uidx) >= 0;
            xu(uidx) <= vehicle.tmax;
            %restrict attitude
%             
%             xu(xidxp1(1)) >= -60*pi/180;
%             xu(xidxp1(1)) <= 60*pi/180;
        end
%         x(N) = A*x(N) + B*u(N-1) + dist;
%         a.k.a x(N) = x(N-1)
          xu(xidxp1) == xu(xidx) ;
    cvx_end
    x_u = reshape(xu,2+n,N);
    vehicle.control_cvx.u = x_u(3:end,:);
    vehicle.control_cvx.x = x_u(1:2,:);
    vehicle.control_cvx.Uff = vehicle.control_cvx.u ;
    vehicle.control_cvx.F = C* vehicle.control_cvx.x + D * vehicle.control_cvx.u;
    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.solved = true;
    vehicle.control_cvx.Ftarget = Ftarget_in;
end

vehicle.control_cvx.iter = vehicle.control_cvx.iter + 1;
if(vehicle.control_cvx.iter < size(vehicle.control_cvx.Uff,2))
    vehicle.U = vehicle.control_cvx.Uff(:,vehicle.control_cvx.iter);
end

end
function vehicle = control_cvx(vehicle, Ftarget)


if((~vehicle.control_cvx.solved) || any(Ftarget ~= vehicle.control_cvx.Ftarget))
    
    %x = Ax * Bu where x = [theta, theta_dot]
    %y = [Fx; delta_Fz] = [-W*theta ; ]
    %y = [-W , 0] * x + [-1 -1 -1 -1] * u
    N = 300;
    n = 4;
    Ftrim = [0 ; -vehicle.weight];
    utrim =  vehicle.weight / 4;
    
    Q = diag([1 10]);
    r = 0.1;
    R = diag(r*ones(n,1));
    
    C = [-vehicle.weight , 0 ; 0, 0];
    D = [0 0 0 0 ; -1 -1 -1 -1];
    
    A = [C'*Q*C , -C'*Q*D ; -D'*Q*C , D'*Q*D + R];
    fd = Ftarget - Ftrim;
    c = -2*[C'*Q; D'*Q] * fd;
    
    Abar = repmat({A},1,N);
    Abar = blkdiag(Abar{:});
    
    cbart = repmat(c,N,1)';
    
    
    cvx_begin
    cvx_quiet(false)
    variable xu(2*N + n*N);
    
    minimize transpose(xu)*Abar*xu + cbart*xu;
    subject to 
        xu(1) == vehicle.x(1);
        xu(2) == vehicle.x(2);
        for i = 1:N-1
            xi = (i-1)*(2+n)+1;
            xidx = xi:(xi+1);
            xidxp1 = xidx + (2+n);
            ui = xi+2; 
            uidx = ui:(ui+n-1);
            %X(k+1) = A*x + B*u
            xu(xidxp1) == vehicle.sysd.a * xu(xidx) + vehicle.sysd.b * xu(uidx);
            %max thrust
            xu(uidx) >= -utrim;
            xu(uidx) <= vehicle.tmax-utrim;
            %restrict attitude
            xu(xidxp1(1)) >= -15*pi/180;
            xu(xidxp1(1)) <= 15*pi/180;
        end
    tic;
    
    cvx_end
   toc;
    x_u = reshape(xu,2+n,N);
    Utrim =utrim * ones(n,1);
    vehicle.control_cvx.u = x_u(3:end,:);
    vehicle.control_cvx.x = x_u(1:2,:);
    vehicle.control_cvx.Uff = vehicle.control_cvx.u + repmat(Utrim,1,N);
    vehicle.control_cvx.F = C* vehicle.control_cvx.x + D * vehicle.control_cvx.u + repmat(Ftrim,1,N);
    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.solved = true;
    vehicle.control_cvx.Ftarget = Ftarget;
end

vehicle.control_cvx.iter = vehicle.control_cvx.iter + 1;
vehicle.U = vehicle.control_cvx.Uff(:,vehicle.control_cvx.iter);

end
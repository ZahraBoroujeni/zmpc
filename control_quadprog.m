function vehicle = control_quadprog(vehicle, Ftarget_in)

%x = Ax * Bu where x = [theta, theta_dot]
%y = [Fx; delta_Fz] = [-W*theta ; ]
%y = [-W , 0] * x + [-1 -1 -1 -1] * u
T = 100;
n = size(vehicle.sysd.a,1);
m = size(vehicle.sysd.b,2);


if(~isfield(vehicle.control_cvx,'x0'));
    vehicle.control_cvx.x0  = zeros(n,1);
    vehicle.control_cvx.x0(2) = -vehicle.weight/vehicle.sysd.c(2,2);
end

x0  = vehicle.control_cvx.x0;
x0(vehicle.sysdthetaIndex) = vehicle.theta;
x0(vehicle.sysdqIndex) = vehicle.q;
x0(vehicle.sysduIndex) = vehicle.u;



Ftarget = Ftarget_in;
theta_Fx = asind(Ftarget(1)/vehicle.weight);
theta_Fx = max(min(theta_Fx,15),-15);
Ftarget(1) = sind(theta_Fx)*vehicle.weight;


Qy = diag([5 ; 50]);
Qyfinal = diag([5; 50]);
r =  1;
R = diag(r*ones(m,1));
C = vehicle.sysd.C;
C(2,:) = C(2,:) / cos(vehicle.theta);
D = vehicle.sysd.D;

    
if(~isfield(vehicle.control_cvx,'H'))
   %build the H matrix 

    A = vehicle.sysd.a ;
    B = vehicle.sysd.b ;
    
    
    Qx = C'*Qy*C ;
    Qxu = -C'*Qy*D;
    Qu = D'*Qy*D + R;
    
    Qxf = C'*Qyfinal*C;

    
    QSR = repmat({[Qx Qxu ; Qxu' Qu]},1,T-1);
    
    vehicle.control_cvx.H = blkdiag(Qu,QSR{:},Qxf);
    vehicle.control_cvx.Qxu = Qxu;
   
     vehicle.control_cvx.fastC = zeros((T+1)*(n),(T)*(n+m));
    I = eye(n);
    vehicle.control_cvx.fastC(1:n,1:(n+m)) = [-B,  I];
    fastC = [-A , -B, eye(n)];
    
    rows = 1:n; rows = rows + n;
    cols = (m+1):(m+1+(2*n+m)-1);
    for i = 2:T
        vehicle.control_cvx.fastC(rows,cols) = fastC;
        rows = rows + n;
        cols = cols + (n+m);
    end
    vehicle.control_cvx.fastC(rows,cols-n-m) = [0*A , -B, (eye(n)-A)];


    
    Fu = [eye(m) ; -eye(m)];
    vehicle.control_cvx.fu = repmat([vehicle.tmax * ones(m,1); vehicle.tmin * ones(m,1)],T,1);
    l = 2*m;
    vehicle.control_cvx.fastP = zeros((T-1)*l,T*(n+m));

    rows = 1:l;
    cols = 1:m;
    for i = 1:(T)
        vehicle.control_cvx.fastP(rows, cols) = Fu;
        rows = rows + l;
        cols = cols + (n+m);
    end
    
    vehicle.control_cvx.problem.H = (vehicle.control_cvx.H+vehicle.control_cvx.H');
    vehicle.control_cvx.problem.Aeq = vehicle.control_cvx.fastC;
    vehicle.control_cvx.problem.lb = -Inf(T*(n+m),1);
    vehicle.control_cvx.problem.ub = Inf(T*(n+m),1);
    for i = 1:T
        uIndex = (1:m) + (i-1)*(n+m);
        vehicle.control_cvx.problem.lb(uIndex) = vehicle.tmin * ones(m,1);
        vehicle.control_cvx.problem.ub(uIndex) = vehicle.tmax * ones(m,1);
    end
    
    vehicle.control_cvx.problem.solver = 'quadprog';
    vehicle.control_cvx.problem.options = optimset('Algorithm','interior-point-convex',...
            'LargeScale','on','Display','off');
    vehicle.control_cvx.problem.options.TolCon = 1e-3;
    vehicle.control_cvx.problem.options.TolFun = 1e-3;
    vehicle.control_cvx.problem.options.TolX = 1e-3;
%problem.options = optimset('Algorithm','trust-region-reflective','Display','off');
% 
%               Algorithm: [ active-set | interior-point | interior-point-convex | levenberg-marquardt | ...
%                            sqp | trust-region-dogleg | trust-region-reflective ]

    vehicle.control_cvx.usol = zeros(m*T,1);
    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.numSol = 0;
    vehicle.control_cvx.deltaT = 0;

end


if( 1 || vehicle.control_cvx.iter == 2 || (~vehicle.control_cvx.solved))

    tstart = tic;
    w = zeros(n,1);
    wdelta = vehicle.sysdMy.b *  vehicle.estimator_dist.Myd;
    w(vehicle.sysdthetaIndex) = wdelta(1);
    w(vehicle.sysdqIndex) = wdelta(2);
    
    b = [vehicle.sysd.A*x0 + w; repmat(w, T-1,1); w];
    
    
    ydesired = Ftarget; %q = 0 [Fx,Fz] = Ftarget
    cx = -2*C'*Qy * ydesired;
    cxf = -2*C'*Qyfinal* ydesired;
    cu = -2*D'*Qy * ydesired;
    
    q = cx; qf = cxf;
    r = cu;
    
    g = [r+2*vehicle.control_cvx.Qxu'*x0 ; q ; repmat([r ; q],T-2,1) ; r ; qf];
    
    %    X = quadprog(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
    %     structure with matrix 'H' in PROBLEM.H, the vector 'f' in PROBLEM.f,
    %     the linear inequality constraints in PROBLEM.Aineq and PROBLEM.bineq,
    %     the linear equality constraints in PROBLEM.Aeq and PROBLEM.beq, the
    %     lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub, the start
    %     point in PROBLEM.x0
    
    vehicle.control_cvx.problem.f = g;
    vehicle.control_cvx.problem.x0 = vehicle.control_cvx.usol;
    vehicle.control_cvx.problem.Beq = b;
    
    [z, ~, ~]  = quadprog(vehicle.control_cvx.problem);
    vehicle.control_cvx.deltaT = vehicle.control_cvx.deltaT + toc(tstart);
    vehicle.control_cvx.numSol = vehicle.control_cvx.numSol + 1;

    u_x = reshape(z,n+m,T);
    usol = u_x(1:m,:);
    xsol = u_x(m+1:end,:);
    
    
    vehicle.control_cvx.usol = usol;
    vehicle.control_cvx.xsol = xsol;
    vehicle.control_cvx.F = C* vehicle.control_cvx.xsol + D * vehicle.control_cvx.usol;
    vehicle.control_cvx.Ftarget = Ftarget_in;

    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.solved = true;
end


vehicle.control_cvx.iter = vehicle.control_cvx.iter + 1;
if(vehicle.control_cvx.iter < size(vehicle.control_cvx.usol,2))
    vehicle.U = vehicle.control_cvx.usol(:,vehicle.control_cvx.iter);
    vehicle.control_cvx.x0 = vehicle.sysd.a *  x0  +  vehicle.sysd.b * vehicle.U;
end

end
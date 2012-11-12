function vehicle = control_cvx_noineq(vehicle, Ftarget_in)

%x = Ax * Bu where x = [theta, theta_dot]
%y = [Fx; delta_Fz] = [-W*theta ; ]
%y = [-W , 0] * x + [-1 -1 -1 -1] * u
N = 50;
n = 4;
x0  = [vehicle.theta; vehicle.q ; vehicle.actD.tstates];

Ftarget = Ftarget_in;
theta_Fx = asind(Ftarget(1)/vehicle.weight);
theta_Fx = max(min(theta_Fx,15),-15);
Ftarget(1) = sind(theta_Fx)*vehicle.weight;


qq = 0;
Qy = diag([qq; 1 ; 50]);
Qyfinal = diag([qq; 1; 100]);
r = 2;
R = diag(r*ones(n,1));
nStates = size(vehicle.sysd.a,1) + size(vehicle.actD.sysd.a,1);
C = zeros(3,nStates);
C(1,2) = 1; %q = q
C(2,1) = -vehicle.weight; %Fx = -theta*mg
D = [zeros(1,4); 
    zeros(1,4);
     -cos(vehicle.theta)*ones(1,n);]; %Fz = -sum(ucmd) -> should be moved as part of C since thrust is now a state!?

delta = 1;



    
if(~isfield(vehicle.control_cvx,'H'))
   %build the H matrix 

    Atheta = vehicle.sysd.a ;
    Btheta = vehicle.sysd.b ;
    Au = vehicle.actD.sysd.a;
    Bu = vehicle.actD.sysd.b;
    Cu = vehicle.actD.sysd.c;
    
    A = [Atheta , Btheta*Cu ;
        zeros(size(Au,1),size(Atheta,2)) , Au];
    B = [zeros(size(Atheta,1),size(Bu,2)) ; Bu];
    
    
    vehicle.control_cvx.c1 = [repmat(zeros(size(A)),1,N-2) , -eye(size(A)) , eye(size(A))] ;
    
    vehicle.control_cvx.BMy = [vehicle.sysdMy.b ; zeros(size(Au,1),1)];
    
    Bbar = repmat({B},1,N-delta);
    Bbar = blkdiag(Bbar{:});
    Z1 = zeros(size(B,1)*delta,size(Bbar,2));
    Z3 = zeros(size(Bbar,1),size(B,2)*delta);
    Z2 = zeros(size(Z1,1),size(Z3,2));
    Bbar = [Z1 , Z2 ; Bbar, Z3];
   
    

    
    Abar = repmat({-A},1,N-1);
    Abar = blkdiag(Abar{:});
    Abar = blkdiag(zeros(size(A)),Abar,zeros(size(A)));
    Abar = Abar(1:end-size(A,1),size(A,2)+1:end);
    Abar = eye(N*size(A)) + Abar;

    
    Abarinv = pinv(Abar);
    Fbar = Abarinv*Bbar;
    
    Qxx = C'*Qy*C ;
    Qxu = C'*Qy*D;
    Quu = D'*Qy*D + R;
    
    Pxx = C'*Qyfinal*C;
    Pxu = Qxu;
    Puu = D'*Qy*D + R;
    
    Qxbar = repmat({Qxx},1,N-1);
    Qxbar = blkdiag(Qxbar{:},Pxx);
    
    Qubar = repmat({Quu},1,N-1);
    Qubar = blkdiag(Qubar{:},Puu);
    
    Qxubar = repmat({Qxu},1,N-1);
    Qxubar = blkdiag(Qxubar{:},Pxu);
   
    
    vehicle.control_cvx.H = Qubar + Fbar'*(Qxbar*Fbar - 2*Qxubar);
    vehicle.control_cvx.Qubar = Qubar;
    vehicle.control_cvx.Qxbar = Qxbar;
    vehicle.control_cvx.Qxubar = Qxubar;
    vehicle.control_cvx.Abarinv = Abarinv;
    vehicle.control_cvx.Fbar = Fbar;
    
    vehicle.control_cvx.problem.H = (vehicle.control_cvx.H+vehicle.control_cvx.H');
    vehicle.control_cvx.problem.Aeq = vehicle.control_cvx.c1 * vehicle.control_cvx.Fbar;
    vehicle.control_cvx.problem.lb = zeros(N*n,1);
    vehicle.control_cvx.problem.ub = (vehicle.tmax)*ones(N*n,1);
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

    vehicle.control_cvx.usol = zeros(n*N,1);
    vehicle.control_cvx.iter = 0;
end


if( 0 || vehicle.control_cvx.iter == 2 || (~vehicle.control_cvx.solved) || any(Ftarget_in ~= vehicle.control_cvx.Ftarget))

    
bMyd = vehicle.control_cvx.BMy *  vehicle.estimator_dist.Myd;
bbar = [x0 ; repmat(bMyd, N-1,1)];


ydesired = [0; Ftarget]; %q = 0 [Fx,Fz] = Ftarget
cx = -2*C'*Qy * ydesired;
cu = -2*D'*Qy * ydesired;
    
cxbart = repmat(cx,N,1)';
cubart = repmat(cu,N,1)';

dbart = (vehicle.control_cvx.Abarinv*bbar)';

CUbart = 2*dbart*(vehicle.control_cvx.Qxbar * vehicle.control_cvx.Fbar - vehicle.control_cvx.Qxubar) + ...
         cxbart*vehicle.control_cvx.Fbar + cubart;
        

     
%    X = quadprog(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%     structure with matrix 'H' in PROBLEM.H, the vector 'f' in PROBLEM.f,
%     the linear inequality constraints in PROBLEM.Aineq and PROBLEM.bineq,
%     the linear equality constraints in PROBLEM.Aeq and PROBLEM.beq, the
%     lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub, the start
%     point in PROBLEM.x0

vehicle.control_cvx.problem.f = CUbart';
vehicle.control_cvx.problem.x0 = vehicle.control_cvx.usol;
vehicle.control_cvx.problem.Beq = -vehicle.control_cvx.c1 * dbart';

u = quadprog(vehicle.control_cvx.problem);


    u(end-3:end) = u(end-7:end-4);
    vehicle.control_cvx.usol = u;
    vehicle.control_cvx.u = reshape(u,n,N);
    
    vehicle.control_cvx.Uff = vehicle.control_cvx.u;
    x = vehicle.control_cvx.Fbar * u + dbart';
    vehicle.control_cvx.x = reshape(x,nStates,N);
    vehicle.control_cvx.F = C* vehicle.control_cvx.x + D * vehicle.control_cvx.u;
    vehicle.control_cvx.Ftarget = Ftarget_in;

    vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.solved = true;
end


vehicle.control_cvx.iter = vehicle.control_cvx.iter + 1;
if(vehicle.control_cvx.iter < size(vehicle.control_cvx.Uff,2))
    vehicle.U = vehicle.control_cvx.Uff(:,vehicle.control_cvx.iter);
end

end
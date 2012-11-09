function vehicle = control_cvx_noineq(vehicle, Ftarget_in)

%x = Ax * Bu where x = [theta, theta_dot]
%y = [Fx; delta_Fz] = [-W*theta ; ]
%y = [-W , 0] * x + [-1 -1 -1 -1] * u
N = 50;
n = 4;

Ftarget = Ftarget_in;
theta_Fx = asind(Ftarget(1)/vehicle.weight);
theta_Fx = max(min(theta_Fx,15),-15);
Ftarget(1) = sind(theta_Fx)*vehicle.weight;


Qx = diag([1 100]);
r = 2;
R = diag(r*ones(n,1));
C = [-vehicle.weight , 0 ; 0, 0];
D = [0 0 0 0 ; -1 -1 -1 -1];

delta = 1;

if(~isfield(vehicle.control_cvx,'H'))
   %build the H matrix 

    A = vehicle.sysd.a ;
    B = vehicle.sysd.b ;
    
    vehicle.control_cvx.c1 = [repmat(zeros(size(A)),1,N-2) , -eye(size(A)) , eye(size(A))] ;
    
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
    
    Qxx = C'*Qx*C ;
    Qxu = C'*Qx*D;
    Quu = D'*Qx*D + R;
    
    Pxx = Qxx;
    Pxu = Qxu;
    Puu = D'*Qx*D + R;
    
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
    vehicle.control_cvx.problem.options = optimset('Algorithm','interior-point-convex','Display','off');    
%problem.options = optimset('Algorithm','trust-region-reflective','Display','off');
% 
%               Algorithm: [ active-set | interior-point | interior-point-convex | levenberg-marquardt | ...
%                            sqp | trust-region-dogleg | trust-region-reflective ]

    vehicle.control_cvx.usol = zeros(n*N,1);
    vehicle.control_cvx.iter = 0;
end


if(vehicle.control_cvx.iter == 2 || (~vehicle.control_cvx.solved) || any(Ftarget_in ~= vehicle.control_cvx.Ftarget))

bMyd = vehicle.sysdMy.b *  vehicle.estimator_dist.Myd;
bbar = [vehicle.x ; repmat(bMyd, N-1,1)];


fd = Ftarget;
cx = -2*C'*Qx * fd;
cu = -2*D'*Qx * fd;
    
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

% % xHx + c'x = norm(Ax - b) = x'A'Ax -2 b'Ax +b'b ; c'=-2b'A -> c=-2*A'*b 
% tic
% A = chol(problem.H);
% b = -inv(A')*problem.f/2; 
% [u, ~] = bvls(A,b,problem.lb,problem.ub);
% toc

% 
%      cvx_begin
%      
%      cvx_quiet(false);
%      variable u(n*N);
%      
%      
%      minimize transpose(u)*problem.H*u/2 + CUbart*u;
%      subject to
%      u >= 0;
%      u <= vehicle.tmax;
%      problem.Aeq * u == problem.Beq;   
%      cvx_end


    u(end-3:end) = u(end-7:end-4);
    vehicle.control_cvx.usol = u;
    vehicle.control_cvx.u = reshape(u,n,N);
    
    vehicle.control_cvx.Uff = vehicle.control_cvx.u;
    x = vehicle.control_cvx.Fbar * u + dbart';
    vehicle.control_cvx.x = reshape(x,2,N);
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
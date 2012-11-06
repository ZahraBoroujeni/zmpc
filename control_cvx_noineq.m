function vehicle = control_cvx_noineq(vehicle, Ftarget_in)

%x = Ax * Bu where x = [theta, theta_dot]
%y = [Fx; delta_Fz] = [-W*theta ; ]
%y = [-W , 0] * x + [-1 -1 -1 -1] * u
N = 300;
n = 4;
Ftrim = [0 ; -vehicle.weight];
utrim =  vehicle.weight / 4;


Ftarget = Ftarget_in;
theta_Fx = asind(Ftarget(1)/vehicle.weight);
theta_Fx = max(min(theta_Fx,15),-15);
Ftarget(1) = sind(theta_Fx)*vehicle.weight;


Qx = diag([1 10]);
r = 0.1;
R = diag(r*ones(n,1));
C = [-vehicle.weight , 0 ; 0, 0];
D = [0 0 0 0 ; -1 -1 -1 -1];


if(~isfield(vehicle.control_cvx,'H'))
   %build the H matrix 

    
    
    A = vehicle.sysd.a ;
    B = vehicle.sysd.b ;
    
    
    
    Bbar = repmat({B},1,N-1);
    Bbar = blkdiag(Bbar{:});
    Bbar = blkdiag(zeros(size(B)),Bbar,zeros(size(B)));
    Bbar = Bbar(1:end-size(B,1),size(B,2)+1:end);
    

    
    Abar = repmat({-A},1,N-1);
    Abar = blkdiag(Abar{:});
    Abar = blkdiag(zeros(size(A)),Abar,zeros(size(A)));
    Abar = Abar(1:end-size(A,1),size(A,2)+1:end);
    Abar = eye(N*size(A)) + Abar;
        
    Abarinv = pinv(Abar);
    Fbar = Abarinv*Bbar;
    
    Qxx = C'*Qx*C ;
    Quu = D'*Qx*D + R;
    Qxu = C'*Qx*D;
    
    Pxx = 10*Qxx;
    Pxu = 10*Qxu;
    Puu = 10*D'*Qx*D + R;
    
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
    
    vehicle.control_cvx.usol = zeros(n*N,1);
    toc;
end


if(0 || (~vehicle.control_cvx.solved) || any(Ftarget_in ~= vehicle.control_cvx.Ftarget))

bMyd = vehicle.sysdMy.b * 1;% vehicle.estimator_dist.Myd;
bbar = [vehicle.x ; repmat(bMyd, N-1,1)];


fd = Ftarget - Ftrim;
cx = -2*C'*Qx * fd;
cu = -2*D'*Qx * fd;
    
cxbart = repmat(cx,N,1)';
cubart = repmat(cu,N,1)';

dbart = (vehicle.control_cvx.Abarinv*bbar)';

CUbart = 2*dbart*(vehicle.control_cvx.Qxbar * vehicle.control_cvx.Fbar - vehicle.control_cvx.Qxubar) + ...
         cxbart*vehicle.control_cvx.Fbar + cubart;
        

    
%      cvx_begin
%      
%      cvx_quiet(false);
%      variable u(n*N);
%      
%      
%      minimize transpose(u)*vehicle.control_cvx.H*u + CUbart*u;
%      subject to
%      u >= -utrim;
%      u <= vehicle.tmax-utrim;
%      cvx_end
 
%    X = quadprog(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%     structure with matrix 'H' in PROBLEM.H, the vector 'f' in PROBLEM.f,
%     the linear inequality constraints in PROBLEM.Aineq and PROBLEM.bineq,
%     the linear equality constraints in PROBLEM.Aeq and PROBLEM.beq, the
%     lower bounds in PROBLEM.lb, the upper bounds in PROBLEM.ub, the start
%     point in PROBLEM.x0
problem.H = (vehicle.control_cvx.H+vehicle.control_cvx.H');
problem.f = CUbart';
problem.lb = -utrim * ones(N*n,1);
problem.ub = (vehicle.tmax-utrim)*ones(N*n,1);
problem.solver = 'quadprog';
problem.options = optimset('Algorithm','interior-point-convex','Display','on');
%problem.options = optimset('Algorithm','trust-region-reflective','Display','off');
% 
%               Algorithm: [ active-set | interior-point | interior-point-convex | levenberg-marquardt | ...
%                            sqp | trust-region-dogleg | trust-region-reflective ]

problem.x0 = vehicle.control_cvx.usol;

% % xHx + c'x = norm(Ax - b) = x'A'Ax -2 b'Ax +b'b ; c'=-2b'A -> c=-2*A'*b 
% tic
% A = chol(problem.H);
% b = -inv(A')*problem.f/2; 
% [u, ~] = bvls(A,b,problem.lb,problem.ub);
% toc

tic;
u = quadprog(problem);
toc;

    u(end-3:end) = u(end-7:end-4);
    vehicle.control_cvx.usol = u;
    Utrim =utrim * ones(n,1);
    vehicle.control_cvx.u = reshape(u,n,N);
    
    vehicle.control_cvx.Uff = vehicle.control_cvx.u + repmat(Utrim,1,N);
    x = vehicle.control_cvx.Fbar * u + dbart';
    vehicle.control_cvx.x = reshape(x,2,N);
    vehicle.control_cvx.F = C* vehicle.control_cvx.x + D * vehicle.control_cvx.u + repmat(Ftrim,1,N);
    vehicle.control_cvx.Ftarget = Ftarget_in;

        vehicle.control_cvx.iter = 0;
    vehicle.control_cvx.solved = true;
end


vehicle.control_cvx.iter = vehicle.control_cvx.iter + 1;
if(vehicle.control_cvx.iter <= size(vehicle.control_cvx.Uff,2))
    vehicle.U = vehicle.control_cvx.Uff(:,vehicle.control_cvx.iter);
end

end
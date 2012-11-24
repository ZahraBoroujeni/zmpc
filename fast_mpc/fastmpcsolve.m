function [U, X] = fastmpcsolve(A,B,C,D,Qy,Qf,R, T, yd, w, x0, U0, X0, ubounds, xbounds)

n = size(B,1);
m = size(B,2);


if(~exist('ubounds','var') || isempty(ubounds))
   ubounds = [-1e6 , 1e6]; 
end

if(~exist('xbounds','var') || isempty(xbounds))
   xbounds = [-1e6 , 1e6]; 
end

if(~exist('X0','var') || isempty(X0))
    X0 = zeros(n,T);
end

if(~exist('U0','var') || isempty(U0))
    U0 = zeros(m,T);
end


Qx = C'*Qy*C;
Qu = R;

r = -2*yd'*Qy*D;
q = -2*yd'*Qy*C;
qf = -2*yd'*Qf*C;


% system description
sys.A = A;
sys.B = B;
sys.xmin =  xbounds(1)*ones(n,1);
sys.xmax =  xbounds(2)*ones(n,1);
sys.umin =  ubounds(1)*ones(m,1);
sys.umax =  ubounds(2)*ones(m,1);
sys.n = n;
sys.m = m;
sys.Q = Qx;
sys.R = Qu;
sys.w = w;
sys.q = q;
sys.qf = qf;
sys.r = r;

% fast MPC parameters
params.T = T; 
params.Qf = Qf;        % final state cost
params.kappa = 0.01;   % barrier parameter
params.niters = 5;     % number of newton steps
params.quiet = false;

[X,U,telapsed] = fmpc_step_old(sys,params,X0,U0,x0);


end


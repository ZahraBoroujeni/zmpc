function [U, X] = fastmpcsolve(A,B,C,D,Qy,Qyf,R, T, yd, w, x0, U0, X0, ubounds, xbounds)

n = size(B,1);
m = size(B,2);


if(~exist('ubounds','var') || isempty(ubounds))
   ubounds = [-1e6 , 1e6]; 
end

if(~exist('xbounds','var') || isempty(xbounds))
   xbounds = [-1e1 , 1e1]; 
end

if(~exist('X0','var') || isempty(X0))
    X0 = zeros(n,T);
end

if(~exist('U0','var') || isempty(U0))
    U0 = zeros(m,T);
end


Qx = C'*Qy*C;
Qxf = C'*Qyf*C;
Qu = R;

r = -2*D'*Qy*yd;
q = -2*C'*Qy*yd;
qf = -2*C'*Qyf*yd;

    
    

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
params.Qf = Qxf;        % final state cost
params.kappa = 0.01;   % barrier parameter
params.niters = 5;     % number of newton steps
params.quiet = true;

%%
%Debugging
% 
% %construct H:
% QSR = repmat({blkdiag(Qx,R)},1,T-1);
% H = blkdiag(R,QSR{:},Qxf);
% 
% nz = numel(X0) + numel(U0);
% 
% z = [];
% for i = 1:size(X0,2)
%     z = [z ; U0(:,i) ; X0(:,i) ];
% end
% 
% g = [repmat([r ; q],T-1,1) ; r ; qf];
% 
% Fu = [eye(m) ; -eye(m)];
% l = 2*m;
% P = zeros((T-1)*l,T*(n+m));
% rows = 1:l;
% cols = 1:m;
% for i = 1:(T)
%     P(rows, cols) = Fu;
%     rows = rows + l;
%     cols = cols + (n+m);
% end
% Pt = P';
% 
% h= repmat([11.92 * ones(m,1); - 0 * ones(m,1)],T,1);
% 
% d  = zeros(size(h));
% for i = 1:numel(h)
%     d(i) = 1/(h(i) -  Pt(:,i)'* z);
% end
% 
% gp = Pt*d;
% gf = 2*H*z + g;
% Hp = P'*diag(d)*diag(d)*P;
% 
% Cmpc = zeros(T*n,(T)*(n+m));
% Cmpc(1:n,1:(n+m)) = [-B,  eye(n)];
% Ctemp = [-A , -B, eye(n)];
% 
% rows = 1:n; rows = rows + n;
% cols = (m+1):(m+1+(2*n+m)-1);
% for i = 2:T
%     Cmpc(rows,cols) = Ctemp;
%     rows = rows + n;
%     cols = cols + (n+m);
% end
% 
% Phi = 2*H + params.kappa * Hp;
% Hess = [ Phi , Cmpc' ;
%     Cmpc , zeros(size(Cmpc,1), size(Cmpc,1))];
% 
% nu = zeros(size(Cmpc,1),1);
% b = [A*x0 ; zeros((T-1)*n,1)];
% rd = gf + params.kappa*gp + Cmpc' * nu;
% rp = Cmpc*z - b;
% 
% dzdnu = -Hess \ [rd ; rp];
% dz = dzdnu(1:nz);
% 
% PhiInv = inv(Phi);
% Y = Cmpc * PhiInv *Cmpc';
% beta = -rp + Cmpc*PhiInv*rd;
% dnu2 = -Y\beta;
% dz2 = PhiInv * (-rd - Cmpc'*dnu2);
%%


[X,U,telapsed] = fmpc_step(sys,params,X0,U0,x0);





end


function K = lqr_finite(A,B,N,Q,R,P)
% 

N = N+1;

error(abcdchk(A,B));
error(abcdchk(Q,B));

if(exist('P','var'))
    error(abcdchk(P,B));
else
    P = Q;
end

[n, m] = size(B);

if( any(size(A) ~= size(Q)) ||...
    any(size(A) ~= size(P)) || ...
    any(size(R) ~= size(zeros(m,m)))) 
    error('incorrect size for Q, P or R');
end

Qbar = repmat({sparse(Q)},1,N);
Qbar = blkdiag(Qbar{:});
Qbar((n*N-n+1):(n*N), (n*N-n+1):(n*N)) = P;

Rbar = repmat({sparse(R)},1,N-1);
Rbar = blkdiag(Rbar{:});

Sx = zeros(n*N,n);
Sx(1:n,:) = eye(n);

for i = 2:N
   Sx(((i-1)*n+1):i*n,:) = A * Sx(((i-2)*n+1):(i-1)*n,:);
end

zIA = zeros(n*N,n);
zIA(n+1:end,1:n) = Sx(1:end-n,:);

Su = zeros(n*N,m*(N-1));
for i = 2:N
    Su(((i-1)*n+1):i*n,1:m) = zIA(((i-1)*n+1):i*n,1:n) * B;
end
for j=2:N-1
   Su((j*n+1):end,((j-1)*m+1):(j*m)) = Su((n+1):(end-(j-1)*n),1:m); 
end
    
Sx = sparse(Sx);
Su = sparse(Su);

H = Su'*Qbar*Su + Rbar;
Ft = Su'*Qbar*Sx;

%using A\b instead of inv(A)
K = H\Ft;
K = full(K(1:m,:));
end
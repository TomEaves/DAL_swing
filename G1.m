function dx = G1(t,x)

% Set parameters
alpha = 0.5;
P = 0.8;
K = 1;

% Set graph
G = graph([1 1 2 3],[2 3 4 4]);
Ps = [-1 1 1 -1]';
Ps = Ps*P;

% Set number of nodes
[N,~] = size(G.Nodes);

A = full(adjacency(G));

% Compute phase differences for use in sin().
phs = zeros(N,N);
for i = 1:N
    for j = 1:N
        phs(i,j) = x(j) - x(i);
    end
end
H = sin(phs);

dx = zeros(2*N,1);

dx(1:N) = x((N+1):2*N);
dx((N+1):2*N) = -alpha*x((N+1):2*N)+Ps+K*diag(A*H');
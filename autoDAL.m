N = 4;          % Number of nodes
alpha = 0.5;    % Nodal Damping
P = 0.8;        % Power scale
K = 1;          % Coupling strength

E0 = 5;             % Current initial energy
dE = 0.1;           % Current gap between energy levels
T = 40;             % Target time
newrun = 1;         % Flag for if this is the first run
prefix = 'dir\0';   % dir = directory to save/load files

Etag = num2str(E0*1e8,'%8.0f'); 
fname = strcat(prefix,Etag,'.mat'); % file name for saved data

if ~newrun                              % continuing a run (or starting a new set of parameters from a previous result)

    load(fname)                         % load saved data

else                                            % starting a new run (no previous data available)

    Ps = [-1 1 1 -1]';                  % set power distribution (vector of +-1)
    G = graph([1 1 2 3],[2 3 4 4]);     % set graph
    steadyphases = [-0.205749941669632; 0.205749941669632; 0.205749941669632; -0.205749941669632]; % set steady phases
    A = adjacency(G);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% set random initial condition %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    y0 = zeros(2*N,1);
    y0(1:N) = steadyphases + randn(N,1);
    y0((N+1):2*N) = randn(N,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% scale random initial condition to have energy E0 %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Dth0s = zeros(N);
    Pot = 0;
    for i = 1:N
        for l = 1:N
            Dth0s(i,l) = steadyphases(l) - steadyphases(i);
            Pot = Pot + 1/4/N*A(i,l)*(y0(l)-y0(i)-Dth0s(i,l))^2;
        end
    end
    Kin = sum(y0((N+1):2*N).^2)/2/N;
    Tot = Kin + Pot;
    y0(1:N) = steadyphases + (y0(1:N)-steadyphases)*sqrt(E0/Tot);
    y0((N+1):2*N) = y0((N+1):2*N)*sqrt(E0/Tot);

    save(fname)

end

countvals = [5 10];

threshold = 3;   % Threshold for E(T) beyond which the dynamics are known to be desynchronised

clear gt

for j = 1:2
    
ct = countvals(j);  % set max number of iterations - 5 and then 10 for 
                    % quick and then more careful anlaysis - allows for 
                    % decisions to be made on dE if all ICs are synchronising

while dE>=1e-5      % keep reducing dE until esp_Ec is met

gt = 0.01; gain = []; resid = [];
count = 0;

while ((gt(end)<threshold) && (count<ct))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% run #countmax iterations (defined in power_grid_DAL) %%%%%%
    %%%%%% or fewer if convergence is reached                   %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [y0t,gtt,~,rt] = power_grid_DAL(E0-dE,T,G,alpha,P*Ps,K,steadyphases,y0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% check final initial condition %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [t,x] = ode45(@uk,[0 T],y0t,options);   
    gainvec = .5*sum(x(:,(N+1):end).^2,2)/N;    % compute gain vector (kinetic energy k(t))

    if gainvec(end)<threshold
        gtt(end) = 0;                           % If synchronised, set k(T)=0 for clearer decision making
    end
    
    y0 = y0t; gt = gtt;                         % Update current 'best' initial condition and current gain k(T)
    gain = [gain gt]; resid = [resid rt];       % Keep track of gains and residuals every #countmax iterations
    count = count + 1;

end     % end loop if deshcyronised condition found, or if count > ct.

if gain(end)>threshold                  % Desynchronisation found   
    E0 = E0 - dE;                       % Lower initial energy
    Etag = num2str(E0*1e8,'%8.0f');     
    fname = strcat(prefix,Etag,'.mat'); % set new file name
    save(fname)                         % save current IC ready for lowering its initial energy and iterating

else                % all ICs synchronised
    load(fname)     % reload initial IC
    clear gt    
    dE = dE/2;      % halve dE
    save(fname)     % resave with new dE and iterate again
end

end % dE < eps_dE

if j == 1    
    dE = 1e-4;      % set very low dE if no desynchronisation ICs found during first run of (#ct * #countmax) iterations
    save(fname)     % save this dE
end

end % both initial and final (more cautious) loop finished
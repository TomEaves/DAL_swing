function [y0,gains,resids] = power_grid_DAL(E0,T,G,alpha,Ps,K,th0s,y0)

% E0: initial energy
% T: time horizon for DAL
% G: graph specification
% alpha: damping parameter
% Ps: power input/output per node
% K: couling strength parameter
% th0s: steady state phases
% y0: initial guess

A = adjacency(G);       % Adjacency matrix
[NN,~] = size(G.Nodes); % Number of nodes

if y0==0                % Use y0=0 to run from random IC
    y0 = zeros(2*NN,1);
    y0(1:NN) = th0s + randn(NN,1);
    y0((NN+1):2*NN) = randn(NN,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% scale initial condition to have energy E0 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dth0s = zeros(NN);
Pot = 0;
for i = 1:NN
    for l = 1:NN
        Dth0s(i,l) = th0s(l) - th0s(i);
        Pot = Pot + 1/4/NN*A(i,l)*(y0(l)-y0(i)-Dth0s(i,l))^2;
    end
end
Kin = sum(y0((NN+1):2*NN).^2)/2/NN;
Tot = Kin + Pot;
y0(1:NN) = th0s + (y0(1:NN)-th0s)*sqrt(E0/Tot);
y0((NN+1):2*NN) = y0((NN+1):2*NN)*sqrt(E0/Tot);

W0 = sum(y0((NN+1):2*NN).^2)/2/NN;  % initial kinetic energy (used to find c)

phs = zeros(NN,NN); % initialise phs
H = phs; M = phs;   % initialise H, M

dt = 0.002;         % timestep for swing equation
epsilon = 0.01;     % initial iteration step-size
tspan = 0:dt:T;     % time vector
NT = length(tspan); % number of time steps
neq = 2*NN;         % total number of equations (2*#nodes)

y = zeros(neq,NT);  % initialise y = [theta_{i=1,N}; omega_{i=1,N}]
Y = zeros(neq,NT);  % initialise Y = [phi_{i=1,N}; eta_{i=1,N}]
F = zeros(neq,4);   % initialise RK4 array

% initialise some algorithm parameters
res = 1;
count = 0;
grad_old = 1;
grad_norm_old = norm(grad_old);

countmax = 20;  % escape after #countmax iterations (e.g. to adjust dE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% perform main iteration %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (abs(res)>1e-8) && (count<countmax)

    count = count + 1;
    
    y(:,1) = y0;    % set IC
    
    loop = 1;
    loopcount = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Use RK4 to integrate the swing equation              %%%%%%
    %%%%%% Fixed step size (NOT ode45) allows for interpolation %%%%%% 
    %%%%%% of theta, omega  onto adjoint variable equation      %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while loop          % continue looping until k(T)_new > k(T)_old
    for i = 2:NT       
        yi = y(:,i-1);
        F(:,1) = forward(yi);                   % 'forward' defined later (swing equation)
        F(:,2) = forward(yi+0.5*dt*F(:,1));        
        F(:,3) = forward(yi+0.5*dt*F(:,2));        
        F(:,4) = forward(yi+dt*F(:,3));        
        y(:,i) = yi + (dt/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4)); 
    end
       
    % compute gain (k(T))

    gainvec = 0.5*sum(y(NN+1:2*NN,:).^2)/NN;    
    loc = length(gainvec);
    gain = gainvec(loc);
        
    if count > 1    % compare against previous gain, if it exists

    if (gain >= 0.9999*gains(count-1)) || (loopcount>49)    

        loop = 0;   % move on to adjoint equations if k(T)_new > k(T)_old

    else                                        % no gain improvement
        
        epsilon = epsilon/1.2;                  % reduce epsilon
        y0 = y0_old;                            % reset IC
        Y0 = Y0_old;                            % reset adjoint at t=0
        c = get_c;                              % find new c
        [y(:,1),~,grad_old,~] = update(y0);     % update IC
        y0 = y(:,1);                            % set IC
        W0 = sum(y0((NN+1):2*NN).^2)/2/NN;      % set k(0)
        loopcount = loopcount + 1;              % stop checking after too many loops

    end

    else % there is no previous gain to compare with
        loop = 0;
    end

    end
    

    gains(count) = gain;            % store gain k(T)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% initialise adjiont variables at t=T %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Y(1:NN,loc) = 0;
    Y(NN+1:2*NN,loc) = y(NN+1:2*NN,loc);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Use RK4 to integrate the adjoin equation from t=T to t=0 %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ode4 adjoint
    dt = -dt;                   % reverse dt (backwards in time)
    for i = loc-1:-1:1      
        yi = y(:,i+1);
        yim1 = y(:,i);
        yimh = 0.5*(yi+yim1);   % interpolate forward variables onto half-step
        Yi = Y(:,i+1);
        F(:,1) = adjoint(yi,Yi);
        F(:,2) = adjoint(yimh,Yi+0.5*dt*F(:,1));        
        F(:,3) = adjoint(yimh,Yi+0.5*dt*F(:,2));        
        F(:,4) = adjoint(yim1,Yi+dt*F(:,3));        
        Y(:,i) = Yi + (dt/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4)); 
    end
    dt = -dt;                   % reset dt
    
    % set ICs
    y0 = y(:,1);
    Y0 = Y(:,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% compute new IC %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    c = get_c;                                      % compute c        
    y0_old = y0;                                    % store IC
    Y0_old = Y0;                                    % store IC
    [y0,res,grad_new,grad_angle] = update(y0_old);  % update IC   
    W0 = sum(y0((NN+1):2*NN).^2)/2/NN;              % set k(0)
    resids(count) = res;                            % store residual             
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% check angle between consecutive gradients %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if count>1  % only perform if previous gradient available

    if (grad_angle > 0.95)                      % closely aligned
        epsilon = 1.25*epsilon;                 % increase step size
        c = get_c;                              % recompute c
        [y0,res,grad_new,~] = update(y0_old);   % recompute IC

    elseif grad_angle < -0.5                    % misaligned
        epsilon = epsilon/2;                    % decrease step size
        c = get_c;                              % recompute c
        [y0,res,grad_new,~] = update(y0_old);   % recompute IC

    end

    end % gradients checked
    
    grad_old = grad_new;                % store gradient
    grad_norm_old = norm(grad_old);     % store gradient size  
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% scale initial condition to have energy E0  %%%%%%
    %%%%% this should already be the case using c    %%%%%%
    %%%%% but sometimes rounding errors when using c %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Pot = 0;
    for i = 1:NN
        for l = 1:NN            
            Pot = Pot + 1/4/NN*A(i,l)*(y0(l)-y0(i)-Dth0s(i,l))^2;
        end
    end    
    Kin = sum(y0((NN+1):2*NN).^2)/2/NN;
    Tot = Kin + Pot;
    y0(1:NN) = th0s + (y0(1:NN)-th0s)*sqrt(E0/Tot);
    y0((NN+1):2*NN) = y0((NN+1):2*NN)*sqrt(E0/Tot);
    W0 = sum(y0((NN+1):2*NN).^2)/2/NN;  % set k(0)
    
end % converged (res<tol) or too many iterations
    



    function dx = forward(x)  
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% the swing equation %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1:NN
            for j = 1:NN
                phs(k,j) = x(j) - x(k);
            end
        end
        H = sin(phs);

        dx = zeros(2*NN,1);
        dx(1:NN) = x((NN+1):2*NN);
        dx((NN+1):2*NN) = -alpha*x((NN+1):2*NN)+Ps+K*diag(A*H');       

    end




    function dX = adjoint(x,X)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% the adjoint equation %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for k = 1:NN
            for j = 1:NN
                phs(k,j) = x(j) - x(k);
            end
        end
        H = cos(phs);

        for j = 1:NN
           M(:,j) = H(:,j) * X(NN+j); 
        end
        dX = zeros(2*NN,1);
        dX(1:NN) = K*X((NN+1):2*NN).*diag(A*H') - K*diag(A*M');
        dX((NN+1):2*NN) = alpha*X((NN+1):2*NN) - X(1:NN);          
        
    end




    function c = get_c
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% compute c %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       disc = -1;
       epsilon = 2*epsilon;
       while disc<0
           epsilon = epsilon/2;

           a1 = 0;
           b1 = 0;
           d1 = 0;     
           S = zeros(NN);
           for k = 1:NN
               for m = 1:NN    
                   for n = 1:NN
                       S(k,m) = S(k,m) + A(m,n)*(y0(n)-y0(m)-Dth0s(m,n))-A(k,n)*(y0(n)-y0(k)-Dth0s(k,n));                        
                   end                    
                   a1 = a1 + A(k,m)*4*epsilon/NN^2*S(k,m)^2;

                   b1 = b1 + A(k,m)*4/NN*(y0(m)-y0(k)-Dth0s(k,m))*S(k,m);
                   b1 = b1 + A(k,m)*4*epsilon/NN^2*(Y0(m)-Y0(k))*S(k,m);                    

                   d1 = d1 + A(k,m)*2/NN*(Y0(m)-Y0(k))*(y0(m)-y0(k)-Dth0s(k,m));
                   d1 = d1 + A(k,m)*epsilon/NN^2*(Y0(m)-Y0(k))^2;                    
               end
           end
            
           a2 = 8*epsilon*W0/NN;
           b2 = -8*W0 - 4*epsilon*Y0((NN+1):2*NN)'*y0((NN+1):2*NN)/NN^2;
           d2 = 2*Y0((NN+1):2*NN)'*y0((NN+1):2*NN)/NN+epsilon*Y0((NN+1):2*NN)'*Y0((NN+1):2*NN)/NN^2;
           
           disc = (b1+b2)^2-4*(a1+a2)*(d1+d2);

       end

       root1 = (-b1-b2 + sqrt((b1+b2)^2-4*(a1+a2)*(d1+d2)))/2/(a1+a2);
       root2 = (-b1-b2 - sqrt((b1+b2)^2-4*(a1+a2)*(d1+d2)))/2/(a1+a2);
       
       % choose the nearest root
       if abs(root1)<=abs(root2)
           c = root1;
       else
           c = root2;
       end
         
    end


    

    function [y0out,res,grad_new,grad_angle] = update(y0in)
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% update IC %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% compute gradients %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        grad_new = Y0/NN;
        for k = 1:NN
            for m = 1:NN
                grad_new(k) = grad_new(k) + 2*c*A(k,m)*(y0in(m)-y0in(k)-Dth0s(k,m))/NN;
            end
        end
        grad_new((NN+1):2*NN,1) = grad_new((NN+1):2*NN,1) - 2*c*y0in((NN+1):2*NN)/NN;
        
        %%%%%%%%%%%%%%%%%%%%
        %%%%%% new IC %%%%%%
        %%%%%%%%%%%%%%%%%%%%

        y0out = y0in + epsilon*grad_new;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% compute redisual %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        res = norm(grad_new./Y0)^2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% compute gradient alignment %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        grad_norm_new = norm(grad_new);
        grad_cross = grad_new'*grad_old;
        grad_angle = grad_cross/grad_norm_new/grad_norm_old;
    
    end

end

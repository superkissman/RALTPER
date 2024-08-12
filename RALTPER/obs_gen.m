function [xV, x_pre_set, xP_pre_set] = obs_gen(inital_pos, dymatic_set, simulation_Time, T, step, warm_time)
% simulation_Time = 10;
% warm_time = 5;
% step = 4;
% T = 0.2;
% inital_pos = [5;5;0;0;0];
% dymatic_set = repmat([0.1;0],[1 (simulation_Time+warm_time)/T]);

n_state=5;      %number of state
n_measurement=2;
q=0.005;    %std of process 
r=0.05;    %std of measurement
Q=q^2*eye(n_state); % covariance of process
R=r^2*eye(n_measurement);        % covariance of measurement  
t = T;
f=@(x)[x(1) + x(4)*cos(x(3))*t;...
       x(2) + x(4)*sin(x(3))*t;...
       x(3) + x(5)*t;...
       x(4);...
       x(5)];  % nonlinear state equations
h=@(x)[x(1);x(2)];                               % measurement equation
s = inital_pos;
x=s+q*randn(n_state,1); %initial state          % initial state with noise
P = eye(n_state);                               % initial state covraiance
N=(simulation_Time+warm_time)/T;                                     % total dynamic steps
xV = zeros(n_state,N);          %estmate        % allocate memory
sV = zeros(n_state,N);          %actual
zV = zeros(n_measurement,N);
%构建真值
for k=1:N
  tmp = [s(1:3);dymatic_set(:,k)];
  s = f(tmp);
  sV(:,k)= s; 
end
xP = zeros(2,2,N)*0.01;
x_pre_set = zeros(n_state,step,N);
xP_pre_set = zeros(2,2,step,N);
%EKF过程
for k=1:N
  s = sV(:,k) + q*randn(n_state,1);
  z = h(s) + r*randn(n_measurement,1);                     % measurments
  zV(:,k)  = z;                             % save measurment
  [x, P] = ekf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  xP(:,:,k)  = P(1:2,1:2);
  x_pre = x;
  xP_pre = P;
  x_pre_set(:,1,k) = x_pre;
  xP_pre_set(:,:,1,k)  = xP_pre(1:2,1:2);
  if(k >= warm_time/T) %给五秒的预跟踪
      for i = 2:step
        [x_pre,xP_pre]=ekf_pre(f,x_pre,xP_pre,Q);
        x_pre_set(:,i,k) = x_pre;
        xP_pre_set(:,:,i,k)  = xP_pre(1:2,1:2);
      end
  end
end

xV = xV(:,warm_time/T+1:end);
x_pre_set = x_pre_set(:,:,warm_time/T+1:end);
xP_pre_set = xP_pre_set(:,:,:,warm_time/T+1:end);

end




function [x1,P]=ekf_pre(fstate,x,P,Q)
    [x1,A]=jaccsd(fstate,x);    %nonlinear update and linearization at current state
    P=A*P*A'+Q; 
end

function [z,A]=jaccsd(fun,x)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)
%
z=fun(x);
n=numel(x);
m=numel(z);
A=zeros(m,n);
h=n*eps;
for k=1:n
    x1=x;
    x1(k)=x1(k)+h*i;
    A(:,k)=imag(fun(x1))/h;
end
end



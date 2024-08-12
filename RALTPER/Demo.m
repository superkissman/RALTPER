clear
close all
clc

% CasADi v3.5.5
addpath('C:\SEU\casadi-windows-matlabR2016a-v3.5.5')

import casadi.*

% 创建 旋转矩阵 && 平移矩阵 匿名函数
Rx = @(x)[cos(x),sin(x);-sin(x),cos(x)]; %顺时针旋转为负
% tx = @(alpha,off,x)[x(1)+cos(alpha)*off(1);x(2)+sin(alpha)*off(2)];
tx = @(x)[x(1);x(2)];

% set the inital and the goal
x0 = [0 ; 1.5 ; 0];   % initial condition.
xs = [28.5 ; 1.5 ; 0]; % Reference posture.


T = 0.2;
sim_tim = 12; % Maximum simulation time
warm_time = 5;
N = 10;
obs_radius = 0.3;

dymatic_set = repmat([1;0],[1 (sim_tim+warm_time)/T]);
[xV, x_pre_set, xP_pre_set] = obs_gen([5;0.4;0;0;0], dymatic_set, sim_tim, T, N, warm_time);

% 定义车身形状
G = [1 0;0 1;-1 0;0 -1];
Width = 3;
Length = 2;
e_g = [Width/2;Length/2;Width/2;Length/2];



v_min = -3.0;v_max = 3.0;
delta_min = -pi/3;delta_max = pi/3;

% state variables
x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

% control variables
v = SX.sym('v'); delta = SX.sym('delta');
controls = [v;delta]; n_controls = length(controls);
% rhs = [v*cos(theta);v*sin(theta);v*tan(delta)/Length]; % system r.h.s
rhs = [v*cos(theta);v*sin(theta);delta]; 

% lagrange variables
lambda = SX.sym('lambda');mu = SX.sym('mu');
lagrange = [lambda;mu];n_lagrange = length(lagrange);

% nonlinear mapping function f(x,u)
f = Function('f',{states,controls},{rhs}); 

% Decision variables (controls)
U = SX.sym('U',n_controls,N-1);

% Decision variables (lambda && mu)
Lambda = SX.sym('Lambda',2,N);
Mu = SX.sym('Mu',4,N);
Nu = SX.sym('Nu',1,N);

% parameters (which include at the initial state of the robot and the reference state)
P = SX.sym('P',n_states, N+1);
A_set = SX.sym('A_set',4,N);
F_set = SX.sym('F_set',4,N);
b_set = SX.sym('b_set',2,N);

% A vector that represents the states over the optimization problem.
X = SX.sym('X',n_states,N);


% Objective function
obj = 0; 

% constraints vector
g = [];  

% weighing matrices (states && controls)
Q = zeros(3,3); Q(1,1) = 0.5;Q(2,2) = 0.5;Q(3,3) = 0.3;
R1 = zeros(2,2); R1(1,1) = 0.5; R1(2,2) = 0.05;
R2 = zeros(2,2); R2(1,1) = 1; R2(2,2) = 0.1;
J = zeros(3,3); J(1,1) = 3;J(2,2) = 3;J(3,3) = 3;


% initial condition constraints
st  = X(:,1); % initial state
g = [g;st-P(:,2)]; 

% the objective function
st_goal = P(:,1);
for k = 1:N-1
    if k==1
        con_last = U(:,k);
    else
        con_last = U(:,k-1);
    end
    st = X(:,k);  con = U(:,k);
    st_next = X(:,k+1);
    st_ref = P(:,k+1); 
    if k <= 4
        % calculate obj 
        obj = obj+(st-st_ref)'*J*(st-st_ref) + con'*R1*con ...
              + (con-con_last)'*R2*(con-con_last) + (st-st_goal)'*Q*(st-st_goal);   
    else
        obj = obj + con'*R1*con + (con-con_last)'*R2*(con-con_last) ...
              + (st-st_goal)'*Q*(st-st_goal);   
    end
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    % add the kinematic constraints
    g = [g;st_next-st_next_euler]; 
end

% Add constraints for collision avoidance
for k = 1:N
    off_m = X(1:2,k);
    alpha_m = X(3,k);
    lambda_m = Lambda(:,k);
    nu_m = Nu(1,k);
    mu_m = Mu(:,k);
    A = [A_set(1:2,k),A_set(3:4,k)];
    F = [F_set(1:2,k),F_set(3:4,k)];
    b = b_set(:,k);
    tmp_o = lambda_m' * A'/(F*2*nu_m);
    % tmp_o = ((F*2*nu_m)\(lambda_m' * A')')';
    tmp1 = nu_m*tmp_o*F*tmp_o';
    tmp2 = lambda_m'*A'*tmp_o';
    g1 = tmp1 - tmp2 -e_g'*mu_m + (tx(off_m)-b)'*lambda_m - nu_m;
    % g1 = -e_g'*mu_m + (tx(off_m)-b)'*lambda_m - nu_m;
    g2 = mu_m'*G + lambda_m'*Rx(alpha_m)';
    g3 = lambda_m'*lambda_m;
    g = [g;g1;g2';g3];
end

% make the decision variable one column  vector
OPT_variables = [reshape(X,3*N,1);reshape(U,2*(N-1),1);reshape(Lambda,2*N,1);reshape(Mu,4*N,1);Nu'];
Params_variables = [reshape(P,3*(N+1),1);reshape(A_set,4*N,1);reshape(F_set,4*N,1);reshape(b_set,2*N,1)];
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', Params_variables);

opts = struct;
opts.ipopt.max_iter = 10000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:3*N) = 0; % equality constraints
args.ubg(1:3*N) = 0; % equality constraints
% g1
args.lbg(3*N+1 : 4 :7*N) = 0.01; % inequality constraints
args.ubg(3*N+1 : 4 :7*N) = inf; % inequality constraints
% g2
args.lbg(3*N+2 : 4 :7*N) = 0; % equality constraints
args.ubg(3*N+2 : 4 :7*N) = 0; % equality constraints
args.lbg(3*N+3 : 4 :7*N) = 0; % equality constraints
args.ubg(3*N+3 : 4 :7*N) = 0; % equality constraints
% g3
args.lbg(3*N+4 : 4 :7*N) = 0; % inequality constraints
args.ubg(3*N+4 : 4 :7*N) = 1; % inequality constraints



args.lbx(1:3:3*N,1) = -10; %state x lower bound
args.ubx(1:3:3*N,1) = 30; %state x upper bound
args.lbx(2:3:3*N,1) = 1; %state y lower bound
args.ubx(2:3:3*N,1) = 2; %state y upper bound
args.lbx(3:3:3*N,1) = -inf; %state theta lower bound
args.ubx(3:3:3*N,1) = inf; %state theta upper bound
args.lbx(3*N+1:2:3*N+2*(N-1),1) = v_min; %v lower bound
args.ubx(3*N+1:2:3*N+2*(N-1),1) = v_max; %v upper bound
args.lbx(3*N+2:2:3*N+2*(N-1),1) = delta_min; %delta lower bound
args.ubx(3*N+2:2:3*N+2*(N-1),1) = delta_max; %delta upper bound
args.lbx(3*N+2*(N-1)+1:1:5*N+2*(N-1),1) = -inf; %lambda && mu lower bound
args.ubx(3*N+2*(N-1)+1:1:5*N+2*(N-1),1) = inf; %lambda && mu upper bound
args.lbx(5*N+2*(N-1)+1:1:10*N+2*(N-1),1) = 0; %lambda && mu && ts lower bound
args.ubx(5*N+2*(N-1)+1:1:10*N+2*(N-1),1) = inf; %lambda && mu && ts upper bound

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;


xx(:,1) = x0; % xx contains the history of states
t(1) = t0;



lambda0 = ones(N,2);
mu0 = ones(N,4);
nu0 = ones(N,1);
u0 = [ones(N-1,1),zeros(N-1,1)];



% Start MPC
mpciter = 1;
record_x = [];
record_u = [];
record_A = [];
record_F = [];
record_b = [];
figure(2)
hold on
axis([-5 30 -1 4])
fill([0,0,30,30],[0,-4,-4,0],'k');
fill([0,30,30,0],[3,3,4,4],'k');

while(norm((x0(1:2)-xs(1:2)),2) > 0.5 && mpciter < sim_tim/T && t0 < sim_tim)
    timer1 = clock;
    
    local_x_pre_set = x_pre_set(:,:,mpciter);
    local_xP_pre_set = xP_pre_set(:,:,:,mpciter);
    [A,F,b] = Creat_AFb(local_x_pre_set, local_xP_pre_set, 0.2, obs_radius, N);
    % A = repmat(A(:,:,1),[1,1,N]);
    % F = repmat(F(:,:,1),[1,1,N]);
    % b = repmat(b(:,1),1,N);
    % Hybrid_A* search the inital path
    % obsInfo = [obs_x,obs_y,obs_radius]
    obsInfo = [xV(1,mpciter);xV(2,mpciter);obs_radius];
    Hybrid_Astar = HyAs1(x0, xs, obsInfo, N);
    timer2 = clock; 
    X0 = Hybrid_Astar;
    if(mpciter == 1)
        P0 = Hybrid_Astar;
    end

    args.p   = [xs;reshape(P0',3*N,1);reshape(A,4*N,1);reshape(F,4*N,1);reshape(b,2*N,1)]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',3*N,1);reshape(u0',2*(N-1),1);reshape(lambda0',2*N,1);reshape(mu0',4*N,1);nu0];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    X0 = reshape(full(sol.x(1:3*N))',3,N)'; % get solution TRAJECTORY
    P0 = [X0(2:end,:);X0(end,:)];
    U_res = reshape(full(sol.x(3*N+1:3*N+2*(N-1)))',2,N-1)';
    lambda0 = reshape(full(sol.x(3*N+2*(N-1)+1:5*N+2*(N-1)))',2,N)';
    mu0 = reshape(full(sol.x(5*N+2*(N-1)+1:9*N+2*(N-1)))',4,N)';
    nu0 = full(sol.x(9*N+2*(N-1)+1:10*N+2*(N-1)));
    [t0, x0, u0] = shift(T, t0, x0, U_res,f);    
    record_x = cat(3,record_x,X0);
    record_u = [record_u;U_res(1,:)];
    record_A = cat(4,record_A,A);
    record_F = cat(4,record_F,F);
    record_b = cat(3,record_b,b);
    timer3 = clock;
    disp(['---The ',num2str(mpciter),'th iteration---']);
    disp(['   The Hybrid_A* time: ',num2str(etime(timer2,timer1))]);
    disp(['   The Optimization time: ',num2str(etime(timer3,timer2))]);
    mpciter = mpciter + 1;
end
Draw_risk_awareness_mul_steps(xs,T,N,G,e_g,record_x,record_u,record_A,record_F,record_b,500);

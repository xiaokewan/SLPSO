%% use the matlab's PSO
% clear
% clc
% fitnessfcn=@slope6;
% nvars=5;            % number of variable
% % para1: air gap, para2: grating period, para3: thickness of silicon
% % grating, para4: thickness of Ag membrane, para5:   wavelength
% lb = [0.5,0.2,0.2,0.02,0.4];
% ub = [2,2,2,0.1,1.6];   
% A = []; b = [];
% Aeq = []; beq = [];
% options = optimoptions('particleswarm','SwarmSize',100,'HybridFcn',@fmincon);
% [x,fval] = particleswarm(fitnessfcn,nvars,A,b,Aeq,beq,lb,ub,options)
% 

%%
% use the standard particle swarm optimization algorithm 


clc;
clear;
close all;

%% Problem Definiton

problem.CostFunction = @slope6 % @(x) Sphere(x);  % Cost Function
problem.nVar = 10;      % Number of Unknown (Decision) Variables
problem.VarMin = -10;   % Lower Bound of Decision Variables
problem.VarMax =  10;   % Upper Bound of Decision Variables

%% Parameters of PSO

% slpso variables
slpso.proSTR = [0.25,0.25,0.25,0.25];			% probabilities of each strategy
slpso.a = 1/6;									% Learning coefficient
slpso.Gs = 10;				 					% learning period
	
% Constriction Coefficients
kappa = 1;
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

params.MaxIt = 1000;        		% Maximum Number of Iterations
params.nSwarm = 5;					% number of swarm
params.nPop = 50;           		% Initial Population Size of each sub-swarm 
params.nPop2 = [40 30 50 49 20];	% Setting Population Size of each sub-swarm
params.w = chi;             		% Intertia Coefficient
params.wdamp = 1;           		% Damping Ratio of Inertia Coefficient
params.c1 = chi*phi1;       		% Personal Acceleration Coefficient
params.c2 = chi*phi2;       		% Social Acceleration Coefficient
params.ShowIterInfo = true; 		% Flag for Showing Iteration Informatin

%% Calling PSO

out = PSO(problem, params);

BestSol = out.BestSol;
BestCosts = out.BestCosts;

%% Results

figure;
% plot(BestCosts, 'LineWidth', 2);
semilogy(BestCosts, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;





		
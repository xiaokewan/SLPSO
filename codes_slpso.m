
%% try to code the SLPSO

function out = PSO(problem, params)

    %% Problem Definiton

    CostFunction = problem.CostFunction;  % Cost Function

    nVar = problem.nVar;        % Number of Unknown (Decision) Variables

    VarSize = [1 nVar];         % Matrix Size of Decision Variables

    VarMin = problem.VarMin;	% Lower Bound of Decision Variables
    VarMax = problem.VarMax;    % Upper Bound of Decision Variables


    %% Parameters of PSO

    MaxIt = params.MaxIt;   % Maximum Number of Iterations

    nPop = params.nPop;     % Population Size (Swarm Size)

    w = params.w;           % Intertia Coefficient
    wdamp = params.wdamp;   % Damping Ratio of Inertia Coefficient
    c1 = params.c1;         % Personal Acceleration Coefficient
    c2 = params.c2;         % Social Acceleration Coefficient
	
	%% Parameters of SLPSO

	proSTR = slpso.proSTR;							% probabilities of each strategy
	a = slpso.a;									% Learning coefficient
	Gs = slpsso.Gs;									% learning period
	S = [];   										% accumulators of 4 strategies
	wei(i) =  log(nPop-i+1)/log( factorial(nPop));	% calculate the weight
	
	% parameters of Roulette wheel selection
	Q = [];			% probalities of strategies
	P = [];			% position in the roulette
	
    % The Flag for Showing Iteration Information
    ShowIterInfo = params.ShowIterInfo;    

    MaxVelocity = 0.2*(VarMax-VarMin);
    MinVelocity = -MaxVelocity;
    
    %% Initialization

    % The Particle Template
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost = inf;

    % Initialize Population Members
    for i=1:nPop

        % Generate Random Solution
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);

        % Update the Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end

    end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt, 1);


    %% Main Loop of PSO
	tgs=1;											% (t%Gs==0)												
	
	for it=1:MaxIt
		while (tgs)
			for i=1:nPop
					
				%% The Roulette wheel sele1 part

				syd = proSTR; 						% probalities of pso stategies
				P=[];
				for j = 1:length(syd)
					P = [P syd(j)/sum(syd)];		% caculate each strategy's probality
				end
				
				Q(1) = P(1);
				
				for o = 2:length(P)
					Q(o) = P(o)+Q(o-1);				% divid the Roulette wheel
				end
				% generate a randm number r
				% select the o if Q（o-1）<=r<Q(o) ture
				flage = 1;							%cycle-index
				count = 1; 	
				z1 = 0;								%strategy 1
				z2 = 0;								%strategy 2 
				z3 = 0;								%strategy 3
				z4 = 0;								%strategy 4

				while(flage)
					r = rand();
					if r<= Q(1)
						z1 = z1+1;
						else if r<= Q(2)
							z2 = z2+1;
						else if r<= Q(3)
							z3 = z3+1;
						else if r<= Q(4) 
							z4 = z4+1;
						end
						end
						end
					end
					
					if count == 100
						flage = 0;
					end
					count = count+1;
				end
				num=[z1 z2 z3 z4];  
				
				[p,index] = max(num);				% get the max num and its position_index
				
				% Update Velocity // select several different stategies
				if index == 1   					% standard pso
				S(1)=S(1)+wei(i);
				
				particle(i).Velocity = w*particle(i).Velocity ...
					+ c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
					+ c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
				end
				
				if index == 2						% fips pso
				S(2)=S(2)+wei(i);
				particles(i).velocity = X*(particles[i].velocity[j] + c * (particles[i].Pm[j] - particles[i].position[j]));
				end
				
				if index == 3						% modified pso only changed the evaluation structure
				S(3)=S(3)+wei(i);
				end
				
				if index == 4						% 
				S(4)=S(4)+wei(i);
				end
				
				% Apply Velocity Limits
				particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
				particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
				
				% Update Position
				particle(i).Position = particle(i).Position + particle(i).Velocity;
				
				% Apply Lower and Upper Bound Limits
				particle(i).Position = max(particle(i).Position, VarMin);
				particle(i).Position = min(particle(i).Position, VarMax);

				% Evaluation
				particle(i).Cost = CostFunction(particle(i).Position);

				% Update Personal Best
				if particle(i).Cost < particle(i).Best.Cost

					particle(i).Best.Position = particle(i).Position;
					particle(i).Best.Cost = particle(i).Cost;

					% Update Global Best
					if particle(i).Best.Cost < GlobalBest.Cost
						GlobalBest = particle(i).Best;
					end            

				end
				% F = [];
				% Evaluate the fitness value of each particle
				% F(i) = GlobalBest.Cost;   
				
			end % end for 1:nPop

			% Store the Best Cost Value
			BestCosts(it) = GlobalBest.Cost;

			% Display Iteration Information
			if ShowIterInfo
				disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
			end

			% Damping Inertia Coefficient
			%w = w * wdamp;

		end %end while
		
		% sort the fitness value F in decending order 
		%
		for i = 1:nPop
			% Find out the jth strategy that generates the particle with fiti'
			S(j) = S(j)+wei(i);
		end
		
		if mod(t,Gs)==0
		
			for j = 1:4
				proSTR'(j)=(1-a)proSTR(j)+aS(j)/Gs;
				S(j) = 0;
			end
			for j=1:4
				proSTR(j)=proSTR'(j)/sum(proSTR');
			end
			
			else if 
				tgs = 0;
				end
				
		end 
		
		t = t+1;
		out.pop = particle;
		out.BestSol = GlobalBest;
		out.BestCosts = BestCosts;
	end % end for 1:ItMax	
end


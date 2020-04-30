	%%
	clc;
	clear;
    close all;
	
	syd = proSTR; 						% probalities of pso stategies
	
	P=zeros(length(syd),1);
	
	for j = 1:length(syd)
		P(j) = syd(j)/sum(syd);		% caculate each strategy's probality
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
				z2 =z2+1;
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


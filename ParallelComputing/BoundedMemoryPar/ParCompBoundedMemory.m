
% this is code for parallel computing using Bounded Memory

% below is the initialization of learning, common to all the workers

times = 100;

% define parameter values

ParameterValues = defineParameters(); 

% solve New Keynesian model steady state

[ SteadyStateValuesNK ] = solveNK_SteadyState( ParameterValues );

% calculate parameter values denoted in the learning tex file

[ ParameterValuesLearning ] = calculateLearningParams( SteadyStateValuesNK, ParameterValues );

% create containers to store firms and households learning data

% actual law of motion, that is, how variables evolve according to the
% model. Include all the variables

% in Parallel case we add a matrix of A instead of just one column as before

ActualLawOfMotion = struct('capital',zeros(1,times),'wage',zeros(1,times),'inflation',zeros(1,times), 'interestRate',zeros(1,times),'markup',zeros(1,times),'consumption',zeros(1,times),'labour',zeros(1,times),'capitalReturn',zeros(1,times),'A',zeros(1,times));

% Household_PLM       - households perceived law of motion 
% HouseholdParameters - households regression (learning) parameters

Household_PLM       = struct('capital',zeros(1,times),'wage',zeros(1,times),'inflation',zeros(1,times), 'interestRate',zeros(1,times),'markup',zeros(1,times),'A',zeros(1,times));
HouseholdParameters = struct('capital_Param',zeros(3,times),'wage_Param',zeros(3,times),'inflation_Param',zeros(3,times), 'interestRate_Param',zeros(3,times),'markup_Param',zeros(3,times));

% Firms_PLM       - firms perceived law of motion 
% FirmsParameters - firms regression (learning) parameters

Firms_PLM       = struct('capital',zeros(1,times),'inflation',zeros(1,times),'markup',zeros(1,times),'A',zeros(1,times));
FirmsParameters = struct('capital_Param',zeros(3,times),'inflation_Param',zeros(3,times),'markup_Param',zeros(3,times));

numPARloops          = 1;
numShockRealizations = 1; % for each parallel loop
% create container to store the data

% initialize each worker at the RE equilibrium by loading RE data

load REmatrix_A

load REvariance

storeResults       = {}; % store results of computations
LearningAlgorithms = {}; % store learning algorithms for households and firms
zMatStore   	   = {}; % store updated zMat in each loop

% these containers are used to store D matrices. It is useful if we want to have different learning algorithms

% households

B_Households       = {};
tmpZ_HH            = {};

% firms

B_Firms            = {};
tmpZ_FF            = {};

E_K = {};
E_S = {};

% memory length of households and firms

% Households

memoryLength_HH = 20; 

% firms

memoryLength_FF = 20; 

for j = 1:numPARloops

	for h = 1:numShockRealizations

		storeResults{j}.Households{h}.PLM        = Household_PLM;
		storeResults{j}.Households{h}.Parameters = HouseholdParameters;
		storeResults{j}.Firms{h}.PLM             = Firms_PLM;
		storeResults{j}.Firms{h}.Parameters      = FirmsParameters;
		storeResults{j}.Actual{h}                = ActualLawOfMotion;

		% create shock realization for each worker

		storeResults{j}.Actual{h}.A(1,1) = 0;
		storeResults{j}.Actual{h}.A(1,2) = 0.005; %abs( sqrt(3)*0.02*(rand-0.5) ); % random uniform variable between -1 and 1 with 0.01 std

		for g = 3:times

			storeResults{j}.Actual{h}.A(1,g) = ParameterValues.rho * storeResults{j}.Actual{h}.A(1,g-1);


		end


		% households parameters

		storeResults{j}.Households{h}.Parameters.capital_Param(:,1)      = [ SteadyStateValuesNK.k*(1-REmatrix_A(1,1)) REmatrix_A(1,1) SteadyStateValuesNK.k*REmatrix_A(2,1) ]';
		storeResults{j}.Households{h}.Parameters.wage_Param(:,1)         = [ SteadyStateValuesNK.w*(1-REmatrix_A(1,2)) REmatrix_A(1,2)*SteadyStateValuesNK.w/SteadyStateValuesNK.k SteadyStateValuesNK.w*REmatrix_A(2,2) ]';
		storeResults{j}.Households{h}.Parameters.interestRate_Param(:,1) = [ SteadyStateValuesNK.R*(1-REmatrix_A(1,3)) REmatrix_A(1,3)*SteadyStateValuesNK.R/SteadyStateValuesNK.k SteadyStateValuesNK.R*REmatrix_A(2,3) ]';
		storeResults{j}.Households{h}.Parameters.inflation_Param(:,1)    = [ 1-REmatrix_A(1,4) REmatrix_A(1,4)/SteadyStateValuesNK.k REmatrix_A(2,4) ]';
		storeResults{j}.Households{h}.Parameters.markup_Param(:,1)       = [ SteadyStateValuesNK.X*(1-REmatrix_A(1,5)) REmatrix_A(1,5)*SteadyStateValuesNK.X/SteadyStateValuesNK.k SteadyStateValuesNK.X*REmatrix_A(2,5) ]';

		% firms parameters

		storeResults{j}.Firms{h}.Parameters.capital_Param(:,1)   = [ SteadyStateValuesNK.k*(1-REmatrix_A(1,1)) REmatrix_A(1,1) SteadyStateValuesNK.k*REmatrix_A(2,1) ]';
		storeResults{j}.Firms{h}.Parameters.inflation_Param(:,1) = [ 1-REmatrix_A(1,4) REmatrix_A(1,4)/SteadyStateValuesNK.k REmatrix_A(2,4) ]';
		storeResults{j}.Firms{h}.Parameters.markup_Param(:,1)    = [ SteadyStateValuesNK.X*(1-REmatrix_A(1,5)) REmatrix_A(1,5)*SteadyStateValuesNK.X/SteadyStateValuesNK.k SteadyStateValuesNK.X*REmatrix_A(2,5) ]';

		% initialize Actual Law of Motion of variables at steady state

		storeResults{j}.Actual{h}.capital(1,1)        = SteadyStateValuesNK.k;
		storeResults{j}.Actual{h}.wage(1,1)           = SteadyStateValuesNK.w;
		storeResults{j}.Actual{h}.inflation(1,1)      = 1;
		storeResults{j}.Actual{h}.interestRate(1,1)   = SteadyStateValuesNK.R;
		storeResults{j}.Actual{h}.markup(1,1)         = SteadyStateValuesNK.X;
		storeResults{j}.Actual{h}.consumption(1,1)    = SteadyStateValuesNK.c;
		storeResults{j}.Actual{h}.labour(1,1)         = SteadyStateValuesNK.L;
		storeResults{j}.Actual{h}.capitalReturn(1,1)  = SteadyStateValuesNK.rk;

		% initialize Perceived Law of Motion at steady state values

		% households

		storeResults{j}.Households{h}.PLM.capital(1,1)      = SteadyStateValuesNK.k;
		storeResults{j}.Households{h}.PLM.wage(1,1)         = SteadyStateValuesNK.w;
		storeResults{j}.Households{h}.PLM.interestRate(1,1) = SteadyStateValuesNK.R;
		storeResults{j}.Households{h}.PLM.inflation(1,1)    = 1;
		storeResults{j}.Households{h}.PLM.markup(1,1)       = SteadyStateValuesNK.X;

		% firms

		storeResults{j}.Firms{h}.PLM.capital(1,1)   = SteadyStateValuesNK.k;
		storeResults{j}.Firms{h}.PLM.inflation(1,1) = 1;
		storeResults{j}.Firms{h}.PLM.markup(1,1)    = SteadyStateValuesNK.X;

		% initialize learning algorithms for household and firm

		% household

		LearningAlgorithms{j}.Households{h}.capital   = BoundedMemoryNeuralNetwork( memoryLength_HH,storeResults{j}.Households{h}.Parameters.capital_Param(:,1),SteadyStateValuesNK.k,SteadyStateValuesNK.k );
		LearningAlgorithms{j}.Households{h}.Wage      = BoundedMemoryNeuralNetwork( memoryLength_HH,storeResults{j}.Households{h}.Parameters.wage_Param(:,1),SteadyStateValuesNK.w,SteadyStateValuesNK.k );
		LearningAlgorithms{j}.Households{h}.Inflation = BoundedMemoryNeuralNetwork( memoryLength_HH,storeResults{j}.Households{h}.Parameters.inflation_Param(:,1),1,SteadyStateValuesNK.k );
		LearningAlgorithms{j}.Households{h}.Interest  = BoundedMemoryNeuralNetwork( memoryLength_HH,storeResults{j}.Households{h}.Parameters.interestRate_Param(:,1),SteadyStateValuesNK.R,SteadyStateValuesNK.k );
		LearningAlgorithms{j}.Households{h}.Markup    = BoundedMemoryNeuralNetwork( memoryLength_HH,storeResults{j}.Households{h}.Parameters.markup_Param(:,1),SteadyStateValuesNK.X,SteadyStateValuesNK.k );

		% initialize values for households learning algos

		LearningAlgorithms{j}.Households{h}.capital.curVarTwo   = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Households{h}.capital.curVarOne   = storeResults{j}.Actual{h}.capital(1,1);

		LearningAlgorithms{j}.Households{h}.Wage.curVarOne      = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Households{h}.Wage.curVarTwo      = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Households{h}.Wage.curVarThree    = storeResults{j}.Actual{h}.wage(1,1);

		LearningAlgorithms{j}.Households{h}.Inflation.curVarOne   = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Households{h}.Inflation.curVarTwo   = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Households{h}.Inflation.curVarThree = storeResults{j}.Actual{h}.inflation(1,1);

		LearningAlgorithms{j}.Households{h}.Interest.curVarOne    = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Households{h}.Interest.curVarTwo    = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Households{h}.Interest.curVarThree  = storeResults{j}.Actual{h}.interestRate(1,1);

		LearningAlgorithms{j}.Households{h}.Markup.curVarOne      = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Households{h}.Markup.curVarTwo      = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Households{h}.Markup.curVarThree    = storeResults{j}.Actual{h}.markup(1,1);

		% update intervals for households

		LearningAlgorithms{j}.Households{h}.Wage.UpdateIntervals();
		LearningAlgorithms{j}.Households{h}.Inflation.UpdateIntervals();
		LearningAlgorithms{j}.Households{h}.Interest.UpdateIntervals();
        LearningAlgorithms{j}.Households{h}.Markup.UpdateIntervals();

		% firms

		LearningAlgorithms{j}.Firms{h}.capital   = BoundedMemoryNeuralNetwork( memoryLength_FF,storeResults{j}.Firms{h}.Parameters.capital_Param(:,1),SteadyStateValuesNK.k,SteadyStateValuesNK.k );
		LearningAlgorithms{j}.Firms{h}.Inflation = BoundedMemoryNeuralNetwork( memoryLength_FF,storeResults{j}.Firms{h}.Parameters.inflation_Param(:,1),1,SteadyStateValuesNK.k );
		LearningAlgorithms{j}.Firms{h}.Markup    = BoundedMemoryNeuralNetwork( memoryLength_FF,storeResults{j}.Firms{h}.Parameters.markup_Param(:,1),SteadyStateValuesNK.X,SteadyStateValuesNK.k );

		% initialize values for firms learning algos


		LearningAlgorithms{j}.Firms{h}.capital.curVarOne   = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Firms{h}.capital.curVarTwo   = storeResults{j}.Actual{h}.A(1,1);
		
		LearningAlgorithms{j}.Firms{h}.Inflation.curVarOne = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Firms{h}.Inflation.curVarTwo = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Firms{h}.Inflation.curVarThree = storeResults{j}.Actual{h}.inflation(1,1);

		LearningAlgorithms{j}.Firms{h}.Markup.curVarOne    = storeResults{j}.Actual{h}.capital(1,1);
		LearningAlgorithms{j}.Firms{h}.Markup.curVarTwo    = storeResults{j}.Actual{h}.A(1,1);
		LearningAlgorithms{j}.Firms{h}.Markup.curVarThree  = storeResults{j}.Actual{h}.markup(1,1);

		% update intervals for firms

		LearningAlgorithms{j}.Firms{h}.Inflation.UpdateIntervals();
		LearningAlgorithms{j}.Firms{h}.Markup.UpdateIntervals();


		% create containers to store zMat and D matrices

		zMatStore{j}{h} = {};
        
        % households

        B_Households{j}{h}       = {};
		tmpZ_HH{j}{h}            = {};

		% firms

		B_Firms{j}{h}            = {};
		tmpZ_FF{j}{h}            = {};  
		E_K{j}{h}                = {};
		E_S{j}{h}                = {};


	end

end

% clear some variables

clear k j g h

% open the parpool

parpool(numPARloops);

% the main learning loop 

parfor j = 1:numPARloops
    
    disp('Parallel Loops: ')
    j
    disp('--')

	for h = 1:numShockRealizations
        
        disp('Realization of Shocks: ')
        h
        disp('--');

		for t = 2:times
            
            disp('Time Periods: ')
            t


			% update the state variable for all workers

			storeResults{j}.Actual{h}.capital(1,t) = exp(storeResults{j}.Actual{h}.A(1,t))*storeResults{j}.Actual{h}.capital(1,t-1)^ParameterValues.alpha*storeResults{j}.Actual{h}.labour(1,t-1)^(1-ParameterValues.alpha)-storeResults{j}.Actual{h}.consumption(1,t-1)+(1-ParameterValues.delta)*storeResults{j}.Actual{h}.capital(1,t-1);

			% update capital in learning algorithms

			LearningAlgorithms{j}.Households{h}.capital.curVarThree = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Firms{h}.capital.curVarThree      = storeResults{j}.Actual{h}.capital(1,t);

			LearningAlgorithms{j}.Households{h}.capital.curVarTwo   = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Firms{h}.capital.curVarTwo        = storeResults{j}.Actual{h}.A(1,t);

			% update intervals for households and firms

			LearningAlgorithms{j}.Households{h}.capital.UpdateIntervals();
			LearningAlgorithms{j}.Firms{h}.capital.UpdateIntervals();

			% define zMatrix for households and firms with updated capital

			zMatStore{j}{h} = [ 1 storeResults{j}.Actual{h}.capital(1,t) storeResults{j}.Actual{h}.A(1,t) ]';

			% update D matrix and parameters for both firms and households

			% households

			[ storeResults{j}.Households{h}.Parameters.capital_Param(:,t) ]      = LearningAlgorithms{j}.Households{h}.capital.do_BM_Perceptron();
			[ storeResults{j}.Households{h}.Parameters.wage_Param(:,t) ]         = LearningAlgorithms{j}.Households{h}.Wage.do_BM_Perceptron();
			[ storeResults{j}.Households{h}.Parameters.inflation_Param(:,t) ]    = LearningAlgorithms{j}.Households{h}.Inflation.do_BM_Perceptron();
			[ storeResults{j}.Households{h}.Parameters.interestRate_Param(:,t) ] = LearningAlgorithms{j}.Households{h}.Interest.do_BM_Perceptron();
			[ storeResults{j}.Households{h}.Parameters.markup_Param(:,t) ]       = LearningAlgorithms{j}.Households{h}.Markup.do_BM_Perceptron();

			% firms

			[ storeResults{j}.Firms{h}.Parameters.capital_Param(:,t) ]   = LearningAlgorithms{j}.Firms{h}.capital.do_BM_Perceptron();
			[ storeResults{j}.Firms{h}.Parameters.inflation_Param(:,t) ] = LearningAlgorithms{j}.Firms{h}.Inflation.do_BM_Perceptron();
			[ storeResults{j}.Firms{h}.Parameters.markup_Param(:,t) ]    = LearningAlgorithms{j}.Firms{h}.Markup.do_BM_Perceptron();

			% compute one step ahead forecast / PLM using updated parameters

			% households

			B_Households{j}{h} = [ 1 0 0; storeResults{j}.Households{h}.Parameters.capital_Param(:,t)'; 0 0 ParameterValues.rho ];

			tmpZ_HH{j}{h} = B_Households{j}{h}*zMatStore{j}{h};

			storeResults{j}.Households{h}.PLM.capital(1,t)      = tmpZ_HH{j}{h}(2,1);
			storeResults{j}.Households{h}.PLM.A(1,t)            = tmpZ_HH{j}{h}(3,1);
			storeResults{j}.Households{h}.PLM.wage(1,t)         = storeResults{j}.Households{h}.Parameters.wage_Param(:,t)'*tmpZ_HH{j}{h};      
			storeResults{j}.Households{h}.PLM.interestRate(1,t) = storeResults{j}.Households{h}.Parameters.interestRate_Param(:,t)'*tmpZ_HH{j}{h};      
			storeResults{j}.Households{h}.PLM.inflation(1,t)    = storeResults{j}.Households{h}.Parameters.inflation_Param(:,t)'*tmpZ_HH{j}{h};
			storeResults{j}.Households{h}.PLM.markup(1,t)       = storeResults{j}.Households{h}.Parameters.markup_Param(:,t)'*tmpZ_HH{j}{h};

			% firms

			B_Firms{j}{h} = [1 0 0; storeResults{j}.Firms{h}.Parameters.capital_Param(:,t)'; 0 0 ParameterValues.rho];

			tmpZ_FF{j}{h} = B_Firms{j}{h}*zMatStore{j}{h};
            
			storeResults{j}.Firms{h}.PLM.capital(1,t)   = tmpZ_FF{j}{h}(2,1);
			storeResults{j}.Firms{h}.PLM.A(1,t)         = tmpZ_FF{j}{h}(3,1);
			storeResults{j}.Firms{h}.PLM.inflation(1,t) = storeResults{j}.Firms{h}.Parameters.inflation_Param(:,t)'*tmpZ_FF{j}{h};
			storeResults{j}.Firms{h}.PLM.markup(1,t)    = storeResults{j}.Firms{h}.Parameters.markup_Param(:,t)'*tmpZ_FF{j}{h};

			% solve for Actual Law of Motion (actually, AlM is always the same but we keep this part for consistency)

			[ E_K{j}{h}, E_S{j}{h} ] = Solve_Expecations_PAR( t, SteadyStateValuesNK, ParameterValuesLearning, ParameterValues, storeResults{j},h );


			storeResults{j}.Actual{h}.inflation(1,t)     = SolveInflation2_PAR( ParameterValues, SteadyStateValuesNK, ParameterValuesLearning, t, storeResults{j}, E_K{j}{h}, E_S{j}{h}, h );


			storeResults{j}.Actual{h}.interestRate(1,t)  = SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl*storeResults{j}.Actual{h}.inflation(1,t)-SteadyStateValuesNK.R*(1-ParameterValues.r_R)*ParameterValues.r_Infl+SteadyStateValuesNK.R;
			storeResults{j}.Actual{h}.markup(1,t)        = -ParameterValues.theta*SteadyStateValuesNK.X/(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)*storeResults{j}.Actual{h}.inflation(1,t)+ParameterValues.theta*SteadyStateValuesNK.X...
    			/(1-ParameterValues.theta)/(1-ParameterValues.theta*ParameterValues.beta)+SteadyStateValuesNK.X/(1-ParameterValues.theta*ParameterValues.beta)*E_K{j}{h}+SteadyStateValuesNK.X;
			storeResults{j}.Actual{h}.capitalReturn(1,t) = storeResults{j}.Actual{h}.interestRate(1,t)/storeResults{j}.Actual{h}.inflation(1,t)-1+ParameterValues.delta;
			storeResults{j}.Actual{h}.labour(1,t)        = storeResults{j}.Actual{h}.capital(1,t)*(storeResults{j}.Actual{h}.markup(1,t)*storeResults{j}.Actual{h}.capitalReturn(1,t)/ParameterValues.alpha/exp(storeResults{j}.Actual{h}.A(1,t)))^(1/(1-ParameterValues.alpha));
			storeResults{j}.Actual{h}.wage(1,t)          = (1-ParameterValues.alpha)*exp(storeResults{j}.Actual{h}.A(1,t))*(storeResults{j}.Actual{h}.capital(1,t)/storeResults{j}.Actual{h}.labour(1,t))^ParameterValues.alpha/storeResults{j}.Actual{h}.markup(1,t);
            storeResults{j}.Actual{h}.consumption(1,t)   = storeResults{j}.Actual{h}.wage(1,t)*storeResults{j}.Actual{h}.labour(1,t)^(1-ParameterValues.eta);


			% update learning algos for the next iteration

			% household

			LearningAlgorithms{j}.Households{h}.capital.curVarOne   = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Households{h}.capital.curVarTwo   = storeResults{j}.Actual{h}.A(1,t);

			LearningAlgorithms{j}.Households{h}.Wage.curVarOne      = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Households{h}.Wage.curVarTwo      = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Households{h}.Wage.curVarThree    = storeResults{j}.Actual{h}.wage(1,t);

			LearningAlgorithms{j}.Households{h}.Inflation.curVarOne   = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Households{h}.Inflation.curVarTwo   = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Households{h}.Inflation.curVarThree = storeResults{j}.Actual{h}.inflation(1,t);

			LearningAlgorithms{j}.Households{h}.Interest.curVarOne    = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Households{h}.Interest.curVarTwo    = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Households{h}.Interest.curVarThree  = storeResults{j}.Actual{h}.interestRate(1,t);
			
			LearningAlgorithms{j}.Households{h}.Markup.curVarOne    = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Households{h}.Markup.curVarTwo    = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Households{h}.Markup.curVarThree  = storeResults{j}.Actual{h}.markup(1,t);

			% update intervals for households

			LearningAlgorithms{j}.Households{h}.Wage.UpdateIntervals();
			LearningAlgorithms{j}.Households{h}.Inflation.UpdateIntervals();
			LearningAlgorithms{j}.Households{h}.Interest.UpdateIntervals();
			LearningAlgorithms{j}.Households{h}.Markup.UpdateIntervals();

			% firms

			LearningAlgorithms{j}.Firms{h}.capital.curVarOne = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Firms{h}.capital.curVarTwo = storeResults{j}.Actual{h}.A(1,t);

			LearningAlgorithms{j}.Firms{h}.Inflation.curVarOne   = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Firms{h}.Inflation.curVarTwo   = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Firms{h}.Inflation.curVarThree = storeResults{j}.Actual{h}.inflation(1,t);

			LearningAlgorithms{j}.Firms{h}.Markup.curVarOne      = storeResults{j}.Actual{h}.capital(1,t);
			LearningAlgorithms{j}.Firms{h}.Markup.curVarTwo      = storeResults{j}.Actual{h}.A(1,t);
			LearningAlgorithms{j}.Firms{h}.Markup.curVarThree    = storeResults{j}.Actual{h}.markup(1,t);

			% update intervals for firms

			LearningAlgorithms{j}.Firms{h}.Inflation.UpdateIntervals();
			LearningAlgorithms{j}.Firms{h}.Markup.UpdateIntervals();


		end

	end

end


% kill the parpool
delete(gcp);

% clear redundant variables
%clearvars -except storeResults

% store variables and write some codes to store the data & take the mean
% from MC simulations


% all the variables are stored in storeResults struct
% compute the mean of each variable

% save the data in HDF5 format

%save('MC_Results_NNt_100_w95_3.mat','storeResults','-v7.3')

%clear all;































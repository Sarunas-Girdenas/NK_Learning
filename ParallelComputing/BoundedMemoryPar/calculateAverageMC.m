function [ average ] = calculateAverageMC(structName)

	% Purpose: calculcate the mean of simulations for the given struct
    
    % NOTE: this function uses partial loading functionality, so make sure
    % that '.mat' file is saved using '-v7.3' flag

	% INPUT: 'structName' - struct we want to load, saved from parallel computation
	% OUTPUT: matrix of average of desired variables

	% first, load the variables of interest

	MatFileObject = matfile(structName);
    
    times = 200;
    
    % we know that struct has one variable: 'storeResults'
    
    % initiate struct
    
    average = struct('capital',zeros(1,times),'inflation',zeros(1,times),'interestRate',zeros(1,times),'markup',zeros(1,times),'capitalReturn',zeros(1,times),'labour',zeros(1,times),'wage',zeros(1,times),'consumption',zeros(1,times));
    
    % create structs for all the variables we need
    
    Tempfile          = {};
    
    capitalTemp       = {};
    inflationTemp     = {};
    interestRateTemp  = {};
    markupTemp        = {};
    capitalReturnTemp = {};
    labourTemp        = {};
    wageTemp          = {};
    consumptionTemp   = {};
    
    index = 0;
    
    for j = 1:length(MatFileObject.storeResults) % if use the length(MatFileObject.storeResults), Matlab load the matobject and it takes more time
        
        Temp0 = MatFileObject.storeResults(1,j);
        
        TempFile = Temp0{1};
        
        for k = 1:length(TempFile.Actual{1})
                
            % take the vector of one simulation result
            
                capitalTemp       = TempFile.Actual{k}.capital;
                inflationTemp     = TempFile.Actual{k}.inflation;
                interestRateTemp  = TempFile.Actual{k}.interestRate;
                markupTemp        = TempFile.Actual{k}.markup;
                capitalReturnTemp = TempFile.Actual{k}.capitalReturn;
                labourTemp        = TempFile.Actual{k}.labour;
                wageTemp          = TempFile.Actual{k}.wage;
                consumptionTemp   = TempFile.Actual{k}.consumption;
                
            % update the average of simulation
            
                average.capital      =   (capitalTemp + index .* average.capital) ./ (index + 1);
                average.inflation    =   (inflationTemp + index .* average.inflation) ./ (index + 1);
                average.interestRate =   (interestRateTemp + index .* average.interestRate) ./ (index + 1);
                average.markup       =   (markupTemp + index .* average.markup) ./ (index + 1);
                average.capitalReturn=   (capitalReturnTemp + index .* average.capitalReturn) ./ (index + 1);
                average.labour       =   (labourTemp + index .* average.labour) ./ (index + 1);
                average.wage         =   (wageTemp + index .* average.wage) ./ (index + 1);
                average.consumption  =   (consumptionTemp + index .* average.consumption) ./ (index + 1);
             
                
                index = index + 1;
                
            end
            
        end
        
    end
    
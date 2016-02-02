function [ average ] = calculateMeanMC(structName)

	% Purpose: calculcate the mean of simulations for the given struct
    
    % NOTE: this function uses partial loading functionality, so make sure
    % that '.mat' file is saved using '-v7.3' flag

	% INPUT: 'structName' - struct we want to load, saved from parallel computation
	% OUTPUT: matrix of average of desired variables

	% first, load the variables of interest

	MatFileObject = matfile(structName);
    
    % we know that struct has one variable: 'storeResults'
    
    % create structs for all the variables we need
    
    capitalTemp       = {};
    inflationTemp     = {};
    interestRateTemp  = {};
    markupTemp        = {};
    capitalReturnTemp = {};
    labourTemp        = {};
    wageTemp          = {};
    consumptionTemp   = {};
    
    
    
    for j = 1:length(MatFileObject.storeResults)
        
        for k = 1:1000
            
            if k == 1
                
                TempFile = MatFileObject.storeResults(1,j);
                
                capitalTemp{j}       = TempFile{1}.Actual{k}.capital;
                inflationTemp{j}     = TempFile{1}.Actual{k}.inflation;
                interestRateTemp{j}  = TempFile{1}.Actual{k}.interestRate;
                markupTemp{j}        = TempFile{1}.Actual{k}.markup;
                capitalReturnTemp{j} = TempFile{1}.Actual{k}.capitalReturn;
                labourTemp{j}        = TempFile{1}.Actual{k}.labour;
                wageTemp{j}          = TempFile{1}.Actual{k}.wage;
                consumptionTemp{j}   = TempFile{1}.Actual{k}.consumption;
                
            else
                
                capitalTemp{j}       = capitalTemp{j} + TempFile{1}.Actual{k}.capital;
                inflationTemp{j}     = inflationTemp{j} + TempFile{1}.Actual{k}.inflation;
                interestRateTemp{j}  = interestRateTemp{j} + TempFile{1}.Actual{k}.interestRate;
                markupTemp{j}        = markupTemp{j} + TempFile{1}.Actual{k}.markup;
                capitalReturnTemp{j} = capitalReturnTemp{j} + TempFile{1}.Actual{k}.capitalReturn;
                labourTemp{j}        = labourTemp{j} + TempFile{1}.Actual{k}.labour;
                wageTemp{j}          = wageTemp{j} + TempFile{1}.Actual{k}.wage;
                consumptionTemp{j}   = consumptionTemp{j} + TempFile{1}.Actual{k}.consumption;
                
            end
            
        end
        
    end
    
    % calculate average 
    
    % initiate struct
    
    average = struct('meanCapital',zeros(1,200),'meanInflation',zeros(1,200),'meanInterest',zeros(1,200),'meanMarkup',zeros(1,200),'meanCapRet',zeros(1,200),'meanLabour',zeros(1,200),'meanWage',zeros(1,200),'meanCons',zeros(1,200));
    
    % 12 is the number of parallel loops
    
    for h = 1:12
        
        if h == 1
            
            average.meanCapital   = capitalTemp{h};
            average.meanInflation = inflationTemp{h};
            average.meanInterest  = interestRateTemp{h};
            average.meanMarkup    = markupTemp{h};
            average.meanCapRet    = capitalReturnTemp{h};
            average.meanLabour    = labourTemp{h};
            average.meanWage      = wageTemp{h};
            average.meanCons      = consumptionTemp{h};
            
        else
            
            average.meanCapital   = average.meanCapital + capitalTemp{h};
            average.meanInflation = average.meanInflation + inflationTemp{h};
            average.meanInterest  = average.meanInterest + interestRateTemp{h};
            average.meanMarkup    = average.meanMarkup + markupTemp{h};
            average.meanCapRet    = average.meanCapRet + capitalReturnTemp{h};
            average.meanLabour    = average.meanLabour + labourTemp{h};
            average.meanWage      = average.meanWage + wageTemp{h};
            average.meanCons      = average.meanCons + consumptionTemp{h};
            
        end
        
    end
    
    % finally, divide by 12000
    
            average.meanCapital   = average.meanCapital / 12000;
            average.meanInflation = average.meanInflation / 12000;
            average.meanInterest  = average.meanInterest / 12000;
            average.meanMarkup    = average.meanMarkup / 12000;
            average.meanCapRet    = average.meanCapRet / 12000;
            average.meanLabour    = average.meanLabour / 12000;
            average.meanWage      = average.meanWage / 12000;
            average.meanCons      = average.meanCons / 12000;
            
    
end

	
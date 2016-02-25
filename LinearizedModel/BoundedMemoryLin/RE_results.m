% this file solves the NK model under RE and get the IRFs of variables

% run dynare file (RE) constant gain learning parameter 0.04

dynare NK1.mod

close all

[REmatrix_A,B_ghu,REvariance] = transition_matrix(oo_);

save REmatrix_A REmatrix_A;

save REvariance REvariance;

% run learning file (the shock is the standard deviation of e_A)

MainLearningFile

REresults = struct('capital',zeros(1,times),'wage',zeros(1,times),'inflation',zeros(1,times), 'interestRate',zeros(1,times),'markup',zeros(1,times),'consumption',zeros(1,times),'labour',zeros(1,times));

% save the IRFs from log-deviation form to actual value form

REresults.capital = ((k_e_A+1).*SteadyStateValuesNK.k)';

REresults.wage = ((w_e_A+1).*SteadyStateValuesNK.w)';

REresults.inflation = (Infl_e_A+1)';

REresults.interestRate = ((R_e_A+1).*SteadyStateValuesNK.R)';

REresults.markup = ((X_e_A+1).*SteadyStateValuesNK.X)';

REresults.consumption = ((c_e_A+1).*SteadyStateValuesNK.c)';

REresults.labour = ((L_e_A+1).*SteadyStateValuesNK.L)';

% plot capital

figure1 = figure('Name','Capital');

axes1 = axes('Parent',figure1);

plot1_1 = plot(ActualLawOfMotion.capital(2:end),'b','LineStyle','--','LineWidth',4);

set(plot1_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot1_2 = plot(REresults.capital,'r','LineWidth',4);

set(plot1_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Capital'); title('Capital');

box('off');

set(axes1,'FontSize',18);

legend1 = legend(axes1,'show');

set(legend1,'FontSize',18,'EdgeColor',[1 1 1],'Location','northeast');

line([1 100],[SteadyStateValuesNK.k SteadyStateValuesNK.k],'LineWidth',1.5,'Color','k')

% plot interest rate

figure2 = figure('Name','Interest Rate');

axes2 = axes('Parent',figure2);

plot2_1 = plot(ActualLawOfMotion.interestRate(2:end),'b','LineStyle','--','LineWidth',4);

set(plot2_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot2_2 = plot(REresults.interestRate,'r','LineWidth',4);

set(plot2_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Interest Rate'); title('Interest Rate');

box('off');

set(axes2,'FontSize',18);

legend2 = legend(axes2,'show');

set(legend2,'FontSize',18,'EdgeColor',[1 1 1],'Location','southeast');

line([1 100],[SteadyStateValuesNK.R SteadyStateValuesNK.R],'LineWidth',1.5,'Color','k')

% plot inflation

figure3 = figure('Name','Inflation');

axes3 = axes('Parent',figure3);

plot3_1 = plot(ActualLawOfMotion.inflation(2:end),'b','LineStyle','--','LineWidth',4);

set(plot3_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot3_2 = plot(REresults.inflation,'r','LineWidth',4);

set(plot3_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Inflation'); title('Inflation');

box('off');

set(axes3,'FontSize',18);

legend3 = legend(axes3,'show');

set(legend3,'FontSize',18,'EdgeColor',[1 1 1],'Location','southeast');

line([1 100],[1 1],'LineWidth',1.5,'Color','k')

% plot consumption

figure4 = figure('Name','Consumption');

axes4 = axes('Parent',figure4);

plot4_1 = plot(ActualLawOfMotion.consumption(2:end),'b','LineStyle','--','LineWidth',4);

set(plot4_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot4_2 = plot(REresults.consumption,'r','LineWidth',4);

set(plot4_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Consumption'); title('Consumption');

box('off');

set(axes4,'FontSize',18);

legend4 = legend(axes4,'show');

set(legend4,'FontSize',18,'EdgeColor',[1 1 1],'Location','northeast');

line([1 100],[SteadyStateValuesNK.c SteadyStateValuesNK.c],'LineWidth',1.5,'Color','k')

% plot labour

figure5 = figure('Name','Labour');

axes5 = axes('Parent',figure5);

plot5_1 = plot(ActualLawOfMotion.labour(2:end),'b','LineStyle','--','LineWidth',4);

set(plot5_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot5_2 = plot(REresults.labour,'r','LineWidth',4);

set(plot5_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Labour'); title('Labour');

box('off');

set(axes5,'FontSize',18);

legend5 = legend(axes5,'show');

set(legend5,'FontSize',18,'EdgeColor',[1 1 1],'Location','northeast');

line([1 100],[SteadyStateValuesNK.L SteadyStateValuesNK.L],'LineWidth',1.5,'Color','k')

% plot wage


figure6 = figure('Name','Wage');

axes6 = axes('Parent',figure6);

plot6_1 = plot(ActualLawOfMotion.wage(2:end),'b','LineStyle','--','LineWidth',4);

set(plot6_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot6_2 = plot(REresults.wage,'r','LineWidth',4);

set(plot6_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Wage'); title('Wage');

box('off');

set(axes6,'FontSize',18);

legend6 = legend(axes6,'show');

set(legend6,'FontSize',18,'EdgeColor',[1 1 1],'Location','northeast');

line([1 100],[SteadyStateValuesNK.w SteadyStateValuesNK.w],'LineWidth',1.5,'Color','k')

% Price mark up


figure7 = figure('Name','Price mark up');

axes7 = axes('Parent',figure7);

plot7_1 = plot(ActualLawOfMotion.markup,'b','LineStyle','--','LineWidth',4);

set(plot7_1,'DisplayName','\gamma^H = 0.04, \gamma^R = 0.04','Color',[0 0 1]);

hold on

plot7_2 = plot(REresults.markup,'r','LineWidth',4);

set(plot7_2,'DisplayName','RE','Color',[1 0 0]);

xlabel('t');  ylabel('Mark up'); title('Mark up');

box('off');

set(axes7,'FontSize',18);

legend7 = legend(axes7,'show');

set(legend7,'FontSize',18,'EdgeColor',[1 1 1],'Location','northeast');

line([1 100],[SteadyStateValuesNK.X SteadyStateValuesNK.X],'LineWidth',1.5,'Color','k')


























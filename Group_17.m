%==========================================================================
% EMPIRICAL METHODS IN FINANCE, HEC Lausanne
% Homework Assignment #5
%--------------------------------------------------------------------------
% Authors: Mathieu BRUSATIN, Florent HUMBERT, Joas BERTINOTTI, Yannick CERUTTI
%==========================================================================
addpath(genpath('Kevin Sheppard Toolbox'))
clear all
 %% QUESTION 3
 
%% =======SIMULATIONS=======

EX3.N = 500;    %500 replications of the experiment
EX3.T = 10000;  %We want to have 10'000 realizations of our RV
rng(10000)      %Setting the seed for reproduceability

EX3.normal = sort(normrnd(0,1,[EX3.N,EX3.T]),2);    %Generating N(0,1) samples
EX3.student = sort(stdtdis_rnd([EX3.N,EX3.T],3),2); %Generating t(3) samples with unit variance
EX3.garch = zeros(EX3.N,EX3.T);
for i=1:EX3.N
    EX3.garch(i,:) = tarch_simulate(10000, [0.000005 .05 .945], 1, 0, 1, 'NORMAL');
end
EX3.garch = sort(EX3.garch,2);
%% =======HILL STATISTICS=======

EX3.q25 = 25;
EX3.q100 = 100;
EX3.q500 = 500;

%normal

EX3.hill.data.n25 = (1/(EX3.q25-1))*(sum(log(EX3.normal(:, end-EX3.q25+1:end))...
    -log(EX3.normal(:, end-EX3.q25+1)),2));
EX3.hill.data.n100 = (1/(EX3.q100-1))*(sum(log(EX3.normal(:, end-EX3.q100+1:end))...
    -log(EX3.normal(:, end-EX3.q100+1)),2));
EX3.hill.data.n500 = (1/(EX3.q500-1))*(sum(log(EX3.normal(:, end-EX3.q500+1:end))...
    -log(EX3.normal(:, end-EX3.q500+1)),2));

%student

EX3.hill.data.s25 = (1/(EX3.q25-1))*(sum(log(EX3.student(:, end-EX3.q25+1:end))...
    -log(EX3.student(:, end-EX3.q25+1)),2));
EX3.hill.data.s100 = (1/(EX3.q100-1))*(sum(log(EX3.student(:, end-EX3.q100+1:end))...
    -log(EX3.student(:, end-EX3.q100+1)),2));
EX3.hill.data.s500 = (1/(EX3.q500-1))*(sum(log(EX3.student(:, end-EX3.q500+1:end))...
    -log(EX3.student(:, end-EX3.q500+1)),2));

%garch

EX3.hill.data.g25 = (1/(EX3.q25-1))*(sum(log(EX3.garch(:, end-EX3.q25+1:end))...
    -log(EX3.garch(:, end-EX3.q25+1)),2));
EX3.hill.data.g100 = (1/(EX3.q100-1))*(sum(log(EX3.garch(:, end-EX3.q100+1:end))...
    -log(EX3.garch(:, end-EX3.q100+1)),2));
EX3.hill.data.g500 = (1/(EX3.q500-1))*(sum(log(EX3.garch(:, end-EX3.q500+1:end))...
    -log(EX3.garch(:, end-EX3.q500+1)),2));

%% =======MEAN&STD=======

EX3.hill.normal = [mean(EX3.hill.data.n25), mean(EX3.hill.data.n100), mean(EX3.hill.data.n500);...
    std(EX3.hill.data.n25), std(EX3.hill.data.n100), std(EX3.hill.data.n500)];
EX3.hill.student = [mean(EX3.hill.data.s25), mean(EX3.hill.data.s100), mean(EX3.hill.data.s500);...
    std(EX3.hill.data.s25), std(EX3.hill.data.s100), std(EX3.hill.data.s500)];
EX3.hill.garch = [mean(EX3.hill.data.g25), mean(EX3.hill.data.g100), mean(EX3.hill.data.g500);...
    std(EX3.hill.data.g25), std(EX3.hill.data.g100), std(EX3.hill.data.g500)];

EX3.hill.nexstd = [mean(EX3.hill.data.n25)/5, mean(EX3.hill.data.n100)/10, mean(EX3.hill.data.n500)/sqrt(500)];
EX3.hill.sexstd = [mean(EX3.hill.data.s25)/5, mean(EX3.hill.data.s100)/10, mean(EX3.hill.data.s500)/sqrt(500)];
EX3.hill.gexstd = [mean(EX3.hill.data.g25)/5, mean(EX3.hill.data.g100)/10, mean(EX3.hill.data.g500)/sqrt(500)];
 %% QUESTION 4
 
%% =======DATA PREPARATION=======
[EX4.ndata, EX4.text, EX4.alldata]...
    = xlsread('DATA_HW5.xls');

%We get the 500 largest elements

EX4.returns = (sort(diff(EX4.ndata)...
        ./EX4.ndata(1:end-1)));

%We calculate the threshold

EX4.threshold = EX4.returns(501)*(-1);

%We then calculate negret

EX4.returns = EX4.returns(1:500)*(-1);


%% =======GENERALIZED PARETO=======

[EX4.GPestim,EX4.GPCI] = gpfit(EX4.returns-EX4.threshold);
EX4.Theta = 0.99;

% Quantile 

EX4.quantile = EX4.returns(500) + ((EX4.GPestim(2)/EX4.GPestim(1))...
    *(((length(EX4.ndata(1:end-1))/EX3.N) *(1-EX4.Theta)).^(-EX4.GPestim(1))-1));

% M the amount to be invested 
EX4.M = 10000 ;
EX4.var = EX4.M * EX4.quantile;

%% QUESTION 5

%% ======= i) ESTIMATION =======

%We recompute the unsorted returns
EX5.returns = (diff(EX4.ndata)...
        ./EX4.ndata(1:end-1));
 
%We compute the AR(1) for mu
EX5.AR1 = fitlm(EX5.returns(1:end-1), EX5.returns(2:end));

%We compute the GARCH(1,1) for sigma
[EX5.GARCH.param,~,EX5.GARCH.var]...
    = tarch(EX5.AR1.Residuals.Raw,1,0,1);

%We calculate the z
EX5.z = (EX5.returns(2:end)-EX5.AR1.Fitted)./(EX5.GARCH.var.^0.5);

%We redo the same as in exercise 4
%We get the 500 largest elements

EX5.z = (sort(EX5.z));

%We calculate the threshold

EX5.threshold = EX5.z(501)*(-1);

%We then calculate negz

EX5.z = EX5.z(1:500)*(-1);


%% ======= i) GENERALIZED PARETO=======

[EX5.GPestim,EX5.GPCI] = gpfit(EX5.z-EX5.threshold);
EX5.Theta = 0.99;

% Quantile 

EX5.quantile = EX5.threshold + ((EX5.GPestim(2)/EX5.GPestim(1))...
    *(((length(EX4.ndata(1:end-1))/EX3.N) *(1-EX5.Theta)).^(-EX5.GPestim(1))-1));

% M the amount to be invested 
EX5.M = 10000 ;
EX5.var = EX5.M * EX5.quantile;

%% ======= ii) FORECASTS =======

%The return forecast
EX5.forecast.ret = EX5.AR1.Coefficients.Estimate(2)...
    *EX5.returns;

%The variance forecast
EX5.forecast.var = EX5.GARCH.param(1)...
    +EX5.GARCH.param(2)*(EX5.AR1.Residuals.Raw(end)^2)...
    +EX5.GARCH.param(3)*(EX5.GARCH.var);

EX5.forecast.quantile = EX5.quantile*((EX5.forecast.var).^(0.5))...
    -EX5.forecast.ret(1:end-1);

%% ======= iii) GRAPHS =======

EX4.quantilegraph(1:length(EX5.forecast.var)) = EX4.quantile;

close all
figure1 = figure
plot([EX5.forecast.quantile])
hold on
plot([EX4.quantilegraph])
xlabel('Date')
ylabel('Quantile')
legend('Conditional Quantile','Unconditional Quantile')
EX5.labels = {'1973','1978','1983','1988','1993','1998','2003','2008','2013','2018'};
set(gca, 'XTick', [1,260,520,780,1040,1300,1560,1820,2080,2340], 'XTickLabel', EX5.labels),xlim([1,2400])
saveas(gcf,'graph.png');

EX5.forecast.quantilemean = mean(EX5.forecast.quantile);

%% =======FUNCTIONS=======

function t = stdtdis_rnd(n,df)
    %BASED ON KEVIN SHEPPARD'S UCSD Toolbox
    %ADAPTED BY MATHIEU BRUSATIN
     z = randn(n);
     x = chi2rnd(df,n);
     t = (z*sqrt(df))./sqrt(x);
    t=t./(sqrt(df/(df-2)));
end

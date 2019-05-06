%==========================================================================
% EMPIRICAL METHODS IN FINANCE, HEC Lausanne
% Homework Assignment #4
%--------------------------------------------------------------------------
% Authors: Mathieu BRUSATIN, Florent HUMBERT, Joas BERTINOTTI, Yannick CERUTTI
%==========================================================================
%% Exercise 1
clear all
addpath(genpath('Kevin Sheppard Toolbox'))
%We begin by importing all the data
[data.ndata, data.text, data.alldata]...
    = xlsread('DATA_HW4.xlsx', 'Feuil1');
%We compute simple returns for stocks and bonds
EX1.returns = diff(data.ndata(:,1:2))./data.ndata(1:end-1,1:2);
EX1.returns(:,3) = data.ndata(2:end,3)./5200;
figure1 = figure;
plot(EX1.returns(:,1))
xlabel('Date')
ylabel('Weekly Returns')
saveas(gcf,'retstocks.png');
figure2=figure;
plot(EX1.returns(:,2))
xlabel('Date')
ylabel('Weekly Returns')
saveas(gcf,'retbonds.png');

%% Exercise 2b
EX2.covMatrix = cov(EX1.returns(:,1:2));
EX2.mu = transpose(mean(EX1.returns(:,1:2)));
EX2.riskFree = mean(EX1.returns(:,3));
EX2.lambda2 = 2;
EX2.e = ones(2,1);
EX2.lambda10 = 10;
EX2.alpha2 = inv(EX2.covMatrix*EX2.lambda2)*(EX2.mu - EX2.e*EX2.riskFree);
EX2.alpha10 = inv(EX2.covMatrix*EX2.lambda10)*(EX2.mu - EX2.e*EX2.riskFree);
EX2.minalpha2 = 1-sum(EX2.alpha2);
EX2.minalpha10 = 1-sum(EX2.alpha10);

%% Exercise 3
%a
%KS tests
[EX3.stocks.h,EX3.stocks.p,EX3.stocks.ksstat,EX3.stocks.cv] = ...
    kstest(EX1.returns(:,1));
[EX3.bonds.h,EX3.bonds.p,EX3.bonds.ksstat,EX3.bonds.cv] = ...
    kstest(EX1.returns(:,2));
%LB test
[EX3.LBs.Q, EX3.LBs.pval] = ljungbox((EX1.returns(:,1)-EX1.returns(:,3)),4);
[EX3.LBb.Q, EX3.LBb.pval] = ljungbox((EX1.returns(:,2)-EX1.returns(:,3)),4);

%b
EX3.stocks.AR1 = fitlm(EX1.returns(1:end-1,1), EX1.returns(2:end,1));
EX3.stocks.hac = hac(EX3.stocks.AR1);
EX3.stocks.robTstat(1) = EX3.stocks.AR1.Coefficients.Estimate(1)/(EX3.stocks.hac(1,1)^0.5);
EX3.stocks.robTstat(2) = EX3.stocks.AR1.Coefficients.Estimate(2)/(EX3.stocks.hac(2,2)^0.5);
EX3.bonds.AR1 = fitlm(EX1.returns(1:end-1,2), EX1.returns(2:end,2));
EX3.bonds.hac = hac(EX3.bonds.AR1);
EX3.bonds.robTstat(1) = EX3.bonds.AR1.Coefficients.Estimate(1)/(EX3.bonds.hac(1,1)^0.5);
EX3.bonds.robTstat(2) = EX3.bonds.AR1.Coefficients.Estimate(2)/(EX3.bonds.hac(2,2)^0.5);
% compute robust p-values
EX3.stocks.pValueR = NaN(2,1);
for i = 1 : 2  % account for two tails!
    if EX3.stocks.robTstat(i) >= 0
        EX3.stocks.pValueR(i,1) = 2 * ...
            (1 - tcdf(EX3.stocks.robTstat(i),886-2) );
    else
        EX3.stocks.pValueR(i,1) = 2 * ...
            tcdf(EX3.stocks.robTstat(i),886-2);
    end
end
EX3.bonds.pValueR = NaN(2,1);
for i = 1 : 2  % account for two tails!
    if EX3.bonds.robTstat(i) >= 0
        EX3.bonds.pValueR(i,1) = 2 * ...
            (1 - tcdf(EX3.bonds.robTstat(i),886-2) );
    else
        EX3.bonds.pValueR(i,1) = 2 * ...
            tcdf(EX3.bonds.robTstat(i),886-2);
    end
end
%}
%c

[~,EX3.stocks.LM.pval, EX3.stocks.LM.stat] = ...
    archtest(EX3.stocks.AR1.Residuals.Raw,4);
[~,EX3.bonds.LM.pval, EX3.bonds.LM.stat] = ...
    archtest(EX3.bonds.AR1.Residuals.Raw,4);

%d

[EX3.stocks.GARCH.coeff,~,EX3.stocks.GARCH.var, EX3.stocks.GARCH.cov] = tarch(EX3.stocks.AR1.Residuals.Raw,1,0,1);
EX3.stocks.GARCH.tstat = EX3.stocks.GARCH.coeff./(diag(EX3.stocks.GARCH.cov).^0.5);
[EX3.bonds.GARCH.coeff,~,EX3.bonds.GARCH.var, EX3.bonds.GARCH.cov] = tarch(EX3.bonds.AR1.Residuals.Raw,1,0,1);
EX3.bonds.GARCH.tstat = EX3.bonds.GARCH.coeff./(diag(EX3.bonds.GARCH.cov).^0.5);
%e
EX3.stocks.forecast = ones(52,1);
EX3.bonds.forecast = ones(52,1);

%The first forecast
EX3.stocks.forecast(1) = EX3.stocks.GARCH.coeff(1)...
    +EX3.stocks.GARCH.coeff(2)*(EX3.stocks.AR1.Residuals.Raw(end)^2)...
    +EX3.stocks.GARCH.coeff(3)*(EX3.stocks.GARCH.var(end));

EX3.bonds.forecast(1) = EX3.bonds.GARCH.coeff(1)...
    +EX3.bonds.GARCH.coeff(2)*(EX3.bonds.AR1.Residuals.Raw(end)^2)...
    +EX3.bonds.GARCH.coeff(3)*(EX3.bonds.GARCH.var(end));


%The other steps
for i=1:51
    EX3.stocks.forecast(i+1) = EX3.stocks.GARCH.coeff(1)...
        +(EX3.stocks.GARCH.coeff(2)+EX3.stocks.GARCH.coeff(3))...
        *EX3.stocks.forecast(i);
    EX3.bonds.forecast(i+1) = EX3.bonds.GARCH.coeff(1)...
        +(EX3.bonds.GARCH.coeff(2)+EX3.bonds.GARCH.coeff(3))...
        *EX3.bonds.forecast(i);
end

figure3 = figure;
plot(EX3.stocks.forecast)
xlabel('Date')
ylabel('Variance')
saveas(gcf,'varstock.png');
figure4=figure;
plot(EX3.bonds.forecast)
xlabel('Date')
ylabel('Variance')
saveas(gcf,'varbond.png');
%% Exercise 4
%a

%Expected returns
EX4.stocks.predictedReturns = EX3.stocks.AR1.Fitted;
EX4.bonds.predictedReturns = EX3.bonds.AR1.Fitted;

%Expected variances
EX4.stocks.predictedVariance = EX3.stocks.GARCH.var;
EX4.bonds.predictedVariance = EX3.bonds.GARCH.var;

%Correlation
EX4.corr = corr(EX3.stocks.AR1.Residuals.Raw, EX3.bonds.AR1.Residuals.Raw);
%EX4.corr = corrcoef(EX4.stocks.predictedVariance, EX4.bonds.predictedVariance);

%Covariance through time with constant correlation
EX4.cov = EX4.corr*(EX4.stocks.predictedVariance.^(0.5))...
    .*(EX4.bonds.predictedVariance.^(0.5));

%b
%We prepare the cumulative averages
EX4.mu(1,:) = cumsum(EX4.stocks.predictedReturns)'./(1:length(EX4.stocks.predictedReturns));
EX4.mu(2,:) = cumsum(EX4.bonds.predictedReturns)'./(1:length(EX4.bonds.predictedReturns));
EX4.riskFree = cumsum(EX1.returns(2:end,3))'./(1:length(EX4.bonds.predictedReturns));

%We now calculate the alphas
EX4.alpha2 = zeros(2,length(EX4.stocks.predictedReturns));
EX4.alpha10 = zeros(2,length(EX4.stocks.predictedReturns));
for i=1:length(EX4.stocks.predictedReturns);
    EX4.alpha2(:,i) = inv([EX4.stocks.predictedVariance(i), EX4.cov(i);...
    EX4.cov(i), EX4.bonds.predictedVariance(i)]*EX2.lambda2)*(EX4.mu(:,i) - EX2.e*EX4.riskFree(i));
end
for i=1:length(EX4.stocks.predictedReturns);
    EX4.alpha10(:,i) = inv([EX4.stocks.predictedVariance(i), EX4.cov(i);...
    EX4.cov(i), EX4.bonds.predictedVariance(i)]*EX2.lambda10)*(EX4.mu(:,i) - EX2.e*EX4.riskFree(i));
end
%c
figure5=figure;
EX4.plot2 = [EX4.alpha2' (1-sum(EX4.alpha2))'];
plot(EX4.plot2)
xlabel('Date')
ylabel('Weights')
legend('Equity','Bonds','Cash')
saveas(gcf,'dynweight2.png');
figure6=figure;
EX4.plot10 = [EX4.alpha10' (1-sum(EX4.alpha10))'];
plot(EX4.plot10)
xlabel('Date')
ylabel('Weights')
legend('Equity','Bonds','Cash')
saveas(gcf,'dynweight10.png');
%d

%For 2b approach
EX4.returns2b2 = [EX2.alpha2; (1-sum(EX2.alpha2))]'*EX1.returns(2:end,:)';
EX4.returns2b10 = [EX2.alpha10; (1-sum(EX2.alpha10))]'*EX1.returns(2:end,:)';
%For 4b approach
EX4.returns4b2 = diag([EX4.alpha2; (1-sum(EX4.alpha2))]'*EX1.returns(2:end,:)')';
EX4.returns4b10 = diag([EX4.alpha10; (1-sum(EX4.alpha10))]'*EX1.returns(2:end,:)')';

%Cumulative returns
%EX4.cumreturns2b2 = cumprod(EX4.returns2b2 + 1);
%EX4.cumreturns2b10 = cumprod(EX4.returns2b10 + 1);
%EX4.cumreturns4b2 = cumprod(EX4.returns4b2 + 1);
%EX4.cumreturns4b10 = cumprod(EX4.returns4b10 + 1);

EX4.cumreturns2b2 = cumsum(log(EX4.returns2b2+1));
EX4.cumreturns2b10 = cumsum(log(EX4.returns2b10+1));
EX4.cumreturns4b2 = cumsum(log(EX4.returns4b2+1));
EX4.cumreturns4b10 = cumsum(log(EX4.returns4b10+1));
figure7=figure;
EX4.plotfinal = [EX4.cumreturns2b2' EX4.cumreturns2b10' EX4.cumreturns4b2'...
    EX4.cumreturns4b10'];
plot(EX4.plotfinal)
xlabel('Date')
ylabel('logReturns')
legend('Static2','Static10','Dynamic2','Dynamic10')
saveas(gcf,'finalplot.png');
%e

EX4.TC2 = (abs(EX4.alpha2(1,2:end)-EX4.alpha2(1,1:end-1))+ ...
    abs(EX4.alpha2( 2,2:end)-EX4.alpha2(2,1:end-1)));
EX4.TC10 = (abs(EX4.alpha10(1,2:end)-EX4.alpha10(1,1:end-1))+ ...
    abs(EX4.alpha10(2,2:end)-EX4.alpha10(2,1:end-1)));

EX4.tau2 = (EX4.cumreturns2b2(end)-EX4.cumreturns4b2(end))/sum(EX4.TC2);
EX4.tau10 = (EX4.cumreturns2b10(end)-EX4.cumreturns4b10(end))/sum(EX4.TC10);
    
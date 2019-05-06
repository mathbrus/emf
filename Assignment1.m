clear all
addpath(genpath('Kevin Sheppard Toolbox'))
%% QUESTION 2
%We import the xls file in a way that lets matlab separate char and doubles
[ndata, text, alldata] = xlsread('DATA_HW1.xlsx');

%daily frequency
%simple returns

dailySimpleReturns.Raw = [NaN(5740,1), NaN(5740,1), NaN(5740,1), NaN(5740,1), NaN(5740,1), NaN(5740,1)];
dailySimpleReturns.Raw(1,:) = 0;


for n=1:6
    for i=2:length(dailySimpleReturns.Raw(:,1))
        dailySimpleReturns.Raw(i,n) = ((ndata(i,n) - ndata(i-1,n))/ndata(i-1,n));
    end
end

dailySimpleReturns.Stocks =  dailySimpleReturns.Raw(:,1);
dailySimpleReturns.USTreasury = dailySimpleReturns.Raw(:,2);
dailySimpleReturns.USCorporate = dailySimpleReturns.Raw(:,3);
dailySimpleReturns.USRealEstate = dailySimpleReturns.Raw(:,4);
dailySimpleReturns.USCommodities = dailySimpleReturns.Raw(:,5);
dailySimpleReturns.Dollar = dailySimpleReturns.Raw(:,6);

%log returns

dailyLogReturns.Raw = [NaN(5740,1), NaN(5740,1), NaN(5740,1), NaN(5740,1), NaN(5740,1), NaN(5740,1)];
dailyLogReturns.Raw(1,:) = 0;

for n=1:6
    for i=2:length(dailyLogReturns.Raw(:,1))
        dailyLogReturns.Raw(i,n) = log(ndata(i,n)/ndata(i-1,n));
    end
end

dailyLogReturns.Stocks =  dailyLogReturns.Raw(:,1);
dailyLogReturns.USTreasury = dailyLogReturns.Raw(:,2);
dailyLogReturns.USCorporate = dailyLogReturns.Raw(:,3);
dailyLogReturns.USRealEstate = dailyLogReturns.Raw(:,4);
dailyLogReturns.USCommodities = dailyLogReturns.Raw(:,5);
dailyLogReturns.Dollar = dailyLogReturns.Raw(:,6);

%weekly frequency
%simple returns

weeklySimpleReturns.Raw = [NaN(1148,1), NaN(1148,1), NaN(1148,1), NaN(1148,1), NaN(1148,1), NaN(1148,1)];
weeklySimpleReturns.Raw(1,:) = 0;

for n=1:6
    for i=2:length(weeklySimpleReturns.Raw(:,1))
        weeklySimpleReturns.Raw(i,n) = ((ndata(5*(i-1)+1,n) - ndata(5*(i-2)+1,n))/ndata(5*(i-2)+1,n));
    end
end

weeklySimpleReturns.Stocks =  weeklySimpleReturns.Raw(:,1);
weeklySimpleReturns.USTreasury = weeklySimpleReturns.Raw(:,2);
weeklySimpleReturns.USCorporate = weeklySimpleReturns.Raw(:,3);
weeklySimpleReturns.USRealEstate = weeklySimpleReturns.Raw(:,4);
weeklySimpleReturns.USCommodities = weeklySimpleReturns.Raw(:,5);
weeklySimpleReturns.Dollar = weeklySimpleReturns.Raw(:,6);

%log returns

weeklyLogReturns.Raw = [NaN(1148,1), NaN(1148,1), NaN(1148,1), NaN(1148,1), NaN(1148,1), NaN(1148,1)];
weeklyLogReturns.Raw(1,:) = 0;

for n=1:6
    for i=2:length(weeklyLogReturns.Raw(:,1))
        weeklyLogReturns.Raw(i,n) = log(ndata(5*(i-1)+1,n)/ndata(5*(i-2)+1,n));
    end
end

weeklyLogReturns.Stocks =  weeklyLogReturns.Raw(:,1);
weeklyLogReturns.USTreasury = weeklyLogReturns.Raw(:,2);
weeklyLogReturns.USCorporate = weeklyLogReturns.Raw(:,3);
weeklyLogReturns.USRealEstate = weeklyLogReturns.Raw(:,4);
weeklyLogReturns.USCommodities = weeklyLogReturns.Raw(:,5);
weeklyLogReturns.Dollar = weeklyLogReturns.Raw(:,6);




%Summary statistics

sumStats.Names = {'Stocks','USTreasury','USCorporate','USRealEstate', 'USCommodities', 'Dollar'};
sumStats.StatNames = {'mean', 'variance','skewness','kurtosis','min','max'};

%Daily
%Simple
sumStats.vars = [dailySimpleReturns.Stocks, dailySimpleReturns.USTreasury, dailySimpleReturns.USCorporate, ...
    dailySimpleReturns.USRealEstate, dailySimpleReturns.USCommodities, dailySimpleReturns.Dollar];
sumStats.dailySimple = NaN(6,6);
sumStats.dailySimple(1,:) = mean(sumStats.vars);
sumStats.dailySimple(2,:) = var(sumStats.vars);
sumStats.dailySimple(3,:) = skewness(sumStats.vars);
sumStats.dailySimple(4,:) = kurtosis(sumStats.vars);
sumStats.dailySimple(5,:) = min(sumStats.vars);
sumStats.dailySimple(6,:) = max(sumStats.vars);

sumStats.dailySimpleTable = mat2dataset(sumStats.dailySimple,'VarNames',sumStats.Names,...
    'ObsNames',sumStats.StatNames);

%Log
sumStats.vars = [dailyLogReturns.Stocks, dailyLogReturns.USTreasury, dailyLogReturns.USCorporate, ...
    dailyLogReturns.USRealEstate, dailyLogReturns.USCommodities, dailyLogReturns.Dollar];
sumStats.dailyLog = NaN(6,6);
sumStats.dailyLog(1,:) = mean(sumStats.vars);
sumStats.dailyLog(2,:) = var(sumStats.vars);
sumStats.dailyLog(3,:) = skewness(sumStats.vars);
sumStats.dailyLog(4,:) = kurtosis(sumStats.vars);
sumStats.dailyLog(5,:) = min(sumStats.vars);
sumStats.dailyLog(6,:) = max(sumStats.vars);

sumStats.dailyLogTable = mat2dataset(sumStats.dailyLog,'VarNames',sumStats.Names,...
    'ObsNames',sumStats.StatNames);

%Weekly
%Simple
sumStats.vars = [weeklySimpleReturns.Stocks, weeklySimpleReturns.USTreasury, weeklySimpleReturns.USCorporate, ...
    weeklySimpleReturns.USRealEstate, weeklySimpleReturns.USCommodities, weeklySimpleReturns.Dollar];
sumStats.weeklySimple = NaN(6,6);
sumStats.weeklySimple(1,:) = mean(sumStats.vars);
sumStats.weeklySimple(2,:) = var(sumStats.vars);
sumStats.weeklySimple(3,:) = skewness(sumStats.vars);
sumStats.weeklySimple(4,:) = kurtosis(sumStats.vars);
sumStats.weeklySimple(5,:) = min(sumStats.vars);
sumStats.weeklySimple(6,:) = max(sumStats.vars);

sumStats.weeklySimpleTable = mat2dataset(sumStats.weeklySimple,'VarNames',sumStats.Names,...
    'ObsNames',sumStats.StatNames);

%Log
sumStats.vars = [weeklyLogReturns.Stocks, weeklyLogReturns.USTreasury, weeklyLogReturns.USCorporate, ...
    weeklyLogReturns.USRealEstate, weeklyLogReturns.USCommodities, weeklyLogReturns.Dollar];
sumStats.weeklyLog = NaN(6,6);
sumStats.weeklyLog(1,:) = mean(sumStats.vars);
sumStats.weeklyLog(2,:) = var(sumStats.vars);
sumStats.weeklyLog(3,:) = skewness(sumStats.vars);
sumStats.weeklyLog(4,:) = kurtosis(sumStats.vars);
sumStats.weeklyLog(5,:) = min(sumStats.vars);
sumStats.weeklyLog(6,:) = max(sumStats.vars);

sumStats.weeklyLogTable = mat2dataset(sumStats.weeklyLog,'VarNames',sumStats.Names,...
    'ObsNames',sumStats.StatNames);

%% QUESTION 3

%First thing to interpret the dates
formatIn = 'dd.mm.yyyy';
date = datenum(text(3:end,1), formatIn);

%We can now add them to our log returns
fields.Names = fieldnames(dailyLogReturns);

for i=2:7
    dailyLogReturns.(fields.Names{i}) = [date, dailyLogReturns.(fields.Names{i})];
end

%For weekly returns we need to pick one date out of 5
n = 1:5:5740;

for i=2:7
    weeklyLogReturns.(fields.Names{i}) = [date(n), weeklyLogReturns.(fields.Names{i})];
end

%3a

for i=2:7
    h(i) = histogram(dailyLogReturns.(fields.Names{i})(:,2));
    hold on
    j(i) = histogram(weeklyLogReturns.(fields.Names{i})(:,2))
    title(fields.Names{i})
    hold off
end

%We sort the returns
for i=2:7
    extremeReturns.dailyLow.(fields.Names{i}) = sortrows(dailyLogReturns.(fields.Names{i}), 2);
    extremeReturns.dailyHigh.(fields.Names{i}) = sortrows(dailyLogReturns.(fields.Names{i}), 2, 'descend');
    extremeReturns.weeklyLow.(fields.Names{i}) = sortrows(weeklyLogReturns.(fields.Names{i}), 2);
    extremeReturns.weeklyHigh.(fields.Names{i}) = sortrows(weeklyLogReturns.(fields.Names{i}), 2, 'descend');
end

%We now display the 5 dates 
formatOut = 'dd.mm.yyyy';
extremeReturns.Names = {'DailyMax','DailyMin','WeeklyMax','WeeklyMin'};
extremeReturns.top5 = [extremeReturns.dailyHigh.Stocks(1:5,1), ...
   extremeReturns.dailyLow.Stocks(1:5,1),...
   extremeReturns.weeklyHigh.Stocks(1:5,1), extremeReturns.weeklyLow.Stocks(1:5,1)];
extremeReturns.table = mat2dataset(extremeReturns.top5,'VarNames',extremeReturns.Names);
datestr(extremeReturns.dailyHigh.Stocks(1:5,1), formatIn);
datestr(extremeReturns.dailyLow.Stocks(1:5,1), formatIn);
datestr(extremeReturns.weeklyHigh.Stocks(1:5,1), formatIn);
datestr(extremeReturns.weeklyLow.Stocks(1:5,1), formatIn);
%3b
%Daily
for i=2:7  
    probabilities.dailyHigh.(fields.Names{i}) = normcdf(extremeReturns.dailyHigh.(fields.Names{i})(1:5,2),...
        sumStats.dailyLog(1,(i-1)), sqrt(sumStats.dailyLog(2,(i-1))));
    thresholds.dailyHigh.(fields.Names{i}) = sumStats.dailyLog(1,(i-1))+3*sqrt(sumStats.dailyLog(2,(i-1)));
    thresholds.dailyHighSums(1, (i-1)) =...
        sum(dailyLogReturns.(fields.Names{i})(:,2)>thresholds.dailyHigh.(fields.Names{i}));
    thresholds.dailyHighProba(1, (i-1)) = thresholds.dailyHighSums(1, (i-1))/...
        length(dailyLogReturns.(fields.Names{i})(:,2));
end
for i=2:7
    probabilities.dailyLow.(fields.Names{i}) = normcdf(extremeReturns.dailyLow.(fields.Names{i})(1:5,2),...
        sumStats.dailyLog(1,(i-1)), sqrt(sumStats.dailyLog(2,(i-1))));
    thresholds.dailyLow.(fields.Names{i}) = sumStats.dailyLog(1,(i-1))-3*sqrt(sumStats.dailyLog(2,(i-1)));
    thresholds.dailyLowSums(1, (i-1)) =...
        sum(dailyLogReturns.(fields.Names{i})(:,2)<thresholds.dailyLow.(fields.Names{i}));
    thresholds.dailyLowProba(1, (i-1)) = thresholds.dailyLowSums(1, (i-1))/...
        length(dailyLogReturns.(fields.Names{i})(:,2));
end
%Weekly
for i=2:7
    probabilities.weeklyHigh.(fields.Names{i}) = normcdf(extremeReturns.weeklyHigh.(fields.Names{i})(1:5,2),...
        sumStats.weeklyLog(1,(i-1)), sqrt(sumStats.weeklyLog(2,(i-1))));
    thresholds.weeklyHigh.(fields.Names{i}) = sumStats.weeklyLog(1,(i-1))+3*sqrt(sumStats.weeklyLog(2,(i-1)));
    thresholds.weeklyHighSums(1, (i-1)) =...
        sum(weeklyLogReturns.(fields.Names{i})(:,2)>thresholds.weeklyHigh.(fields.Names{i}));
    thresholds.weeklyHighProba(1, (i-1)) = thresholds.weeklyHighSums(1, (i-1))/...
        length(weeklyLogReturns.(fields.Names{i})(:,2));
end
for i=2:7
    probabilities.weeklyLow.(fields.Names{i}) = normcdf(extremeReturns.weeklyLow.(fields.Names{i})(1:5,2),...
        sumStats.weeklyLog(1,(i-1)), sqrt(sumStats.weeklyLog(2,(i-1))));
    thresholds.weeklyLow.(fields.Names{i}) = sumStats.weeklyLog(1,(i-1))-3*sqrt(sumStats.weeklyLog(2,(i-1)));
    thresholds.weeklyLowSums(1, (i-1)) =...
        sum(weeklyLogReturns.(fields.Names{i})(:,2)<thresholds.weeklyLow.(fields.Names{i}));
    thresholds.weeklyLowProba(1, (i-1)) = thresholds.weeklyLowSums(1, (i-1))/...
        length(weeklyLogReturns.(fields.Names{i})(:,2));
end

%3c
%Daily
warning('off','stats:jbtest:PTooSmall')
for i=2:7
    jarqueBera.Daily.(fields.Names{i}) = [NaN NaN];
    [jarqueBera.Daily.(fields.Names{i})(1), jarqueBera.Daily.(fields.Names{i})(2)] ...
        = jbtest(dailyLogReturns.(fields.Names{i})(:,2),0.05);
end
%Weekly
for i=2:7
    jarqueBera.Weekly.(fields.Names{i}) = [NaN NaN];
    [jarqueBera.Weekly.(fields.Names{i})(1), jarqueBera.Weekly.(fields.Names{i})(2)] ...
        = jbtest(weeklyLogReturns.(fields.Names{i})(:,2),0.05);
end

%3d
%Daily
for i=2:7
    autoCorr.Daily.(fields.Names{i}) = sacf(dailyLogReturns.(fields.Names{i})(:,2), 10, 1);
end
for i=2:7
    [autoCorr.Daily.LB.(fields.Names{i})(:,1), autoCorr.Daily.LB.(fields.Names{i})(:,2)] = ...
        ljungbox(dailyLogReturns.(fields.Names{i})(:,2), 10);
end 

%Weekly
for i=2:7
    autoCorr.Weekly.(fields.Names{i}) = sacf(weeklyLogReturns.(fields.Names{i})(:,2), 10, 1);
end
for i=2:7
    [autoCorr.Weekly.LB.(fields.Names{i})(:,1), autoCorr.Weekly.LB.(fields.Names{i})(:,2)] = ...
        ljungbox(weeklyLogReturns.(fields.Names{i})(:,2), 10);
end
close all
        
%Squared
%Compute returns
for i=2:7
dailyLogReturns.Squared.(fields.Names{i}) = (((dailyLogReturns.(fields.Names{i})(:,2))+1).^2)-1;
weeklyLogReturns.Squared.(fields.Names{i}) = (((weeklyLogReturns.(fields.Names{i})(:,2))+1).^2)-1;
end
for i=2:7
    [autoCorr.Daily.Squared.(fields.Names{i})(:,1), autoCorr.Daily.Squared.(fields.Names{i})(:,2)] = ...
        ljungbox(dailyLogReturns.Squared.(fields.Names{i})(:,1), 10);
end
for i=2:7
    [autoCorr.Weekly.Squared.(fields.Names{i})(:,1), autoCorr.Weekly.Squared.(fields.Names{i})(:,2)] = ...
        ljungbox(weeklyLogReturns.Squared.(fields.Names{i})(:,1), 10);
end

%% QUESTION 4

portfolio.Daily.Returns = zeros(5740,1);
for i=2:7
portfolio.Daily.Returns = portfolio.Daily.Returns + (1/6)*(dailySimpleReturns.(fields.Names{i}));
end

portfolio.Weekly.Returns = zeros(1148,1);
for i=2:7
portfolio.Weekly.Returns = portfolio.Weekly.Returns + (1/6)*(weeklySimpleReturns.(fields.Names{i}));
end

%4a

sumStats.Portfolio = NaN(6,2);
sumStats.Portfolio(1,1) = mean(portfolio.Daily.Returns);
sumStats.Portfolio(2,1) = var(portfolio.Daily.Returns);
sumStats.Portfolio(3,1) = skewness(portfolio.Daily.Returns);
sumStats.Portfolio(4,1) = kurtosis(portfolio.Daily.Returns);
sumStats.Portfolio(5,1) = min(portfolio.Daily.Returns);
sumStats.Portfolio(6,1) = max(portfolio.Daily.Returns);
sumStats.Portfolio(1,2) = mean(portfolio.Weekly.Returns);
sumStats.Portfolio(2,2) = var(portfolio.Weekly.Returns);
sumStats.Portfolio(3,2) = skewness(portfolio.Weekly.Returns);
sumStats.Portfolio(4,2) = kurtosis(portfolio.Weekly.Returns);
sumStats.Portfolio(5,2) = min(portfolio.Weekly.Returns);
sumStats.Portfolio(6,2) = max(portfolio.Weekly.Returns);

sumStats.portfolioTable = mat2dataset(sumStats.Portfolio,'VarNames',{'Daily', 'Weekly'},...
    'ObsNames',sumStats.StatNames);
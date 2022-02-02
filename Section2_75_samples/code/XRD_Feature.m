%% Inport ras data
list = dir('*.ras'); 
% addpath 'D:\Katayama Lab\matlab';


load 'MeasurePara.mat';

% Option 1 for previous data
opts1 = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts1.DataLines = [319, 3319];
opts1.Delimiter = " ";

% Specify column names and types
opts1.VariableNames = ["VarName1", "VarName2"];
opts1.VariableTypes = ["double", "double"];

% Specify file level properties
opts1.ExtraColumnsRule = "ignore";
opts1.EmptyLineRule = "read";
opts1.ConsecutiveDelimitersRule = "join";
opts1.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts1 = setvaropts(opts1, ["VarName1", "VarName2"], "TrimNonNumeric", true);
opts1 = setvaropts(opts1, ["VarName1", "VarName2"], "ThousandsSeparator", ",");

% Option 2 for new data
opts2 = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts2.DataLines = [319, 5319];
opts2.Delimiter = " ";

% Specify column names and types
opts2.VariableNames = ["VarName1", "VarName2"];
opts2.VariableTypes = ["double", "double"];

% Specify file level properties
opts2.ExtraColumnsRule = "ignore";
opts2.EmptyLineRule = "read";
opts2.ConsecutiveDelimitersRule = "join";
opts2.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts2 = setvaropts(opts2, ["VarName1", "VarName2"], "TrimNonNumeric", true);
opts2 = setvaropts(opts2, ["VarName1", "VarName2"], "ThousandsSeparator", ",");

% Apply options for data respectively
figure
for i = 1 : 28
    i
    data1(:, :, i) = readmatrix(list(i).name, opts1);
    plot(data1(:,1,i),data1(:,2,i))
    hold on
end

for i = 29 : length(list)
    i
    data2(:, :, i - 28) = readmatrix(list(i).name, opts2);
    plot(data2(:,1,i-28),data2(:,2,i-28))
    
end
hold off
% Create data matrix
MeasureData1.x = data1(:, 1, 1);
for i = 1 : size(data1, 3)
    MeasureData1.Data(:, i) = data1(:, 2, i);
end

MeasureData2.x = data2(:, 1, 1);
for i = 1 : size(data2, 3)
    MeasureData2.Data(:, i) = data2(:, 2, i);
end

MeasureData1.Name = "XRD";
MeasureData1.xLabel = "Degree";

MeasureData2.Name = "XRD";
MeasureData2.xLabel = "Degree";

%% Create Test Data
id = 26;
if id <= 28
    clear testData
    testData(:, 1) = MeasureData1.x;
    testData(:, 2) = MeasureData1.Data(:, id);
else 
    clear testData
    testData(:, 1) = MeasureData2.x;
    testData(:, 2) = MeasureData2.Data(:, id - 28);
end

figure(1);
scatter(testData(:, 1), testData(:, 2),'filled');

%% Range Selection
degReg = MeasurePara.XRD.degree_region;
degIdx = find((testData(:, 1) >= degReg(1)) & (testData(:, 1) <= degReg(2)));
x = testData(degIdx, 1);
y = testData(degIdx, 2);
figure(2);
scatter(x, y, 'filled');

%% Background Elimination
windowSize = MeasurePara.XRD.WindowSize;
stepSize = MeasurePara.XRD.StepSize;
qtValue = MeasurePara.XRD.QuantileValue;

[backCorrSpec] = msbackadj(x, y, 'WindowSize', windowSize, ...
    'StepSize', stepSize, 'RegressionMethod', 'spline', ...
    'EstimationMethod', 'quantile', 'QuantileValue', ...
    qtValue, 'SmoothMethod', 'loess');

nanIdx = logical(sum(isnan(backCorrSpec),2));
backCorrSpec(nanIdx) = [ ];
x(nanIdx) = [ ];

figure(3);
scatter(x, y, 'filled');
hold on;
scatter(x, backCorrSpec, 'filled');
hold off;
y = backCorrSpec;

%% Smoothing Curve
smthSpan = MeasurePara.XRD.smoothing_span;
yy = smooth(x, y, smthSpan, 'sgolay');

figure(4)
scatter(x, y, 'filled');
hold on;
scatter(x, yy, 'filled');


%% Find Positive Peaks
pks = [ ];
locs = [ ];
peakReg = MeasurePara.XRD.find_peak_region;
peakDis = MeasurePara.XRD.minimum_peak_distance;

for i = 1 : length(peakDis)
    peakRestr = peakReg(i, :);
    degRestr = find((x >= peakRestr(1)) & (x <= peakRestr(2)));
    xT = x(degRestr);
    yyT = yy(degRestr);
    minDis = peakDis(i);
%     [pksT, locsT, wT, ~] = findpeaks(yyT, xT, 'MinPeakDistance', minDis);
[pksT,locsT] = max(yyT);
locsT = xT(locsT);
    pks = [pks; pksT];
    locs = [locs; locsT];
end

scatter(locs, pks, 'filled', '^');
hold off
%% Collecting Features
feat = [ ];
figure(5);hold on; figure(6);hold on
for i = 1 : length(list)
    testData = [ ];
    i
    if i <= 28
        testData(:, 1) = MeasureData1.x;
        testData(:, 2) = MeasureData1.Data(:, i);
    else 
        testData(:, 1) = MeasureData2.x;
        testData(:, 2) = MeasureData2.Data(:, i - 28);
    end
    
    degIdx = find((testData(:, 1) >= degReg(1)) & (testData(:, 1) <= degReg(2)));
    x = testData(degIdx, 1);
    y = testData(degIdx, 2);

    % Eliminate background
    [backCorrSpec] = msbackadj(x, y, 'WindowSize', windowSize, ...
      'StepSize', stepSize, 'RegressionMethod', 'spline', ...
      'EstimationMethod', 'quantile', 'QuantileValue', ...
      qtValue, 'SmoothMethod', 'loess');
    nanIdx = logical(sum(isnan(backCorrSpec),2));
    backCorrSpec(nanIdx) = [];
    x(nanIdx) = [];
    y = backCorrSpec;
figure(5);plot(x,y);xlabel('2\theta [degree]');ylabel('Intensity [-]')
    % Smoothing
    yy = smooth(x, y, smthSpan, 'sgolay');
    figure(6);plot(x,yy);xlabel('2\theta [degree]');ylabel('Intensity [-]')
    % Find Peaks
    pks = [ ];
    locs = [ ];
    pks0 = [ ];
    locs0 = [ ];
    for j = 1 : length(peakDis)
        peakRestr = peakReg(j, :);
        degRestr = find((x >= peakRestr(1)) & (x <= peakRestr(2)));
        xT = x(degRestr);
        yyT = yy(degRestr);
        minDis = peakDis(j);
        %         [pksT, locsT, wT, ~] = findpeaks(yyT, xT, 'MinPeakDistance', minDis);
        [pksT,locsT] = max(yyT);
        locsT = xT(locsT);
        pks0 = [pks0; pksT];
        locs0 = [locs0; locsT];
    end
    
    % Select 11 peaks
    if length(pks0) > 11
            pks0 = [pks0(1); pks0(3 : length(pks0))];
            locs0 = [locs0(1); locs0(3 : length(locs0))];
    end         
        
    pks = [pks; pks0];
    locs = [locs; locs0];
        
    temp = [pks', locs'];    
    feat = [feat; temp];
end
    figure(5);adjfig;figure(6);adjfig
%% Create features table
featData = array2table(feat);
featData.Properties.VariableNames = {'pk1_int', 'pk2_int', 'pk3_int', ...
     'pk4_int', 'pk5_int', 'pk6_int', 'pk7_int', 'pk8_int', 'pk9_int', ...
      'pk1_loc', 'pk2_loc', 'pk3_loc', 'pk4_loc', ...
     'pk5_loc', 'pk6_loc', 'pk7_loc', 'pk8_loc', 'pk9_loc'};
 
 FeatureData_XRD = featData;
for i =1:size(FeatureData_XRD,2)
    FeatureData_XRD.Properties.VariableNames{i} = ['XRD_', FeatureData_XRD.Properties.VariableNames{i}];
end
 
% addpath 'D:\Katayama Lab\matlab';
% load 'Features'; 
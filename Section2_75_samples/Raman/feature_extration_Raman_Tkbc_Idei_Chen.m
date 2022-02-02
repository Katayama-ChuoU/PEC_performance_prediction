% Reading Raman data of Tokubuchi Idei Chen
list = dir('*.CSV');

load('MeasurePara_Raman.mat')
opts = delimitedTextImportOptions("NumVariables", 4);

% 範囲と区切り記号の指定
opts.DataLines = [10, Inf];
opts.Delimiter = " ";

% 列名と型の指定
opts.VariableNames = ["Data", "Var2", "VarName3", "Var4"];
opts.SelectedVariableNames = ["Data", "VarName3"];
opts.VariableTypes = ["double", "string", "double", "string"];

% ファイル レベルのプロパティを指定
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% 変数プロパティを指定
opts = setvaropts(opts, ["Var2", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var4"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["Data", "VarName3"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["Data", "VarName3"], "ThousandsSeparator", ",");

% opts2 = delimitedTextImportOptions("NumVariables", 10);
% 
% % 範囲と区切り記号の指定
% opts2.DataLines = [10, Inf];
% opts2.Delimiter = " ";
% 
% % 列名と型の指定
% opts2.VariableNames = ["Data", "Var2", "VarName3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10"];
% opts2.SelectedVariableNames = ["Data", "VarName3"];
% opts2.VariableTypes = ["double", "string", "double", "string", "string", "string", "string", "string", "string", "string"];
% 
% % ファイル レベルのプロパティを指定
% opts2.ExtraColumnsRule = "ignore";
% opts2.EmptyLineRule = "read";
% opts2.ConsecutiveDelimitersRule = "join";
% opts2.LeadingDelimitersRule = "ignore";
% 
% % 変数プロパティを指定
% opts2 = setvaropts(opts2, ["Var2", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10"], "WhitespaceRule", "preserve");
% opts2 = setvaropts(opts2, ["Var2", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10"], "EmptyFieldRule", "auto");
% opts2 = setvaropts(opts2, ["Data", "VarName3"], "TrimNonNumeric", true);
% opts2 = setvaropts(opts2, ["Data", "VarName3"], "ThousandsSeparator", ",");


for i = 1:length(list)
    i
%     if i <= 28

    Data = readtable(list(i).name,opts); % raw data
    
%     elseif i >= 29
%      Data = readtable(list(i).name,opts2); 
%     end
    data(:,:,i) = table2array(Data);
end
%%



for i=1:size(data,3)
    MeasureData.x(:,i) = data(:,1,i);
    MeasureData.Data(:,i) = data(:,2,i);
end

% remove x column if x(wavelength, 2theta...) have same values in each column
for j = 1:length(list)-1
    x_id(:,j) = MeasureData.x(:,j+1) - MeasureData.x(:,j);
end

if sum(x_id(:)) == 0
    MeasureData.x = MeasureData.x(:,1);
end
MeasureData.Name = "Raman";
MeasureData.xLabel = "Ramanshift";

%% Parameters and conditions check using a single data

    %%
    % tst_dat = readmatrix(list(1).name, opts);
    DataNum = 18;
    tst_dat(:,1) = MeasureData.x; tst_dat(:,2) = MeasureData.Data(:,DataNum);
    figure(1)
    scatter(tst_dat(:,1), tst_dat(:,2), 'filled');

    %% Range selection
    Ramanshift_region = MeasurePara_Raman.Raman_region;
    Removal_region = MeasurePara_Raman.Removal_Ramanshift_region;
    tot_idx = ((tst_dat(:,1) >= Ramanshift_region(1)) & (tst_dat(:,1) <= Ramanshift_region(2)));
%     rem_idx1 = find(((tst_dat(:,1) < Removal_region(1,1)) | (tst_dat(:,1) > Removal_region(1,2))) );
%     rem_idx2 = find(((tst_dat(:,1) < Removal_region(2,1)) | (tst_dat(:,1) > Removal_region(2,2))));
%     [rem_idx, ~] = intersect(rem_idx1,rem_idx2);
%     [tot_idx, ~] = intersect(raman_idx,rem_idx);
    
    figure(2)
    scatter(tst_dat(tot_idx,1), tst_dat(tot_idx,2), 'filled');
    x = tst_dat(tot_idx,1); y = tst_dat(tot_idx,2);
    %% background elmination--------If the data is UV-Vis, skip / comment out this section.-----
    % Background could affect the PEC
    % Some toolbox(Bioinformatics?) is necessary to run this part, install
    % that if necessary.
    
    if MeasureData.Name ~= "UV/Vis"
        
        WindowSize = MeasurePara_Raman.WindowSize;
        StepSize = MeasurePara_Raman.StepSize;
        QuantileValue = MeasurePara_Raman.QuantileValue;
        
        [BackCorrSpec] = msbackadj(x,y,'WindowSize',WindowSize,'StepSize',StepSize,.....
            'RegressionMethod','spline','EstimationMethod', 'quantile', 'QuantileValue', QuantileValue,'SmoothMethod', 'loess');
        
        % remove Nan(not a number)to avoid error
        nan_idx = logical(sum(isnan(BackCorrSpec),2));
        BackCorrSpec(nan_idx) = [];
        x(nan_idx) = [];
        
        figure(3)
        scatter(x,y,'filled')
        hold on
        scatter(x,BackCorrSpec,'filled')
        hold off
        y = BackCorrSpec;
    end
    %% Smoothing curve
    smth_spn = MeasurePara_Raman.smoothing_span;
    yy = smooth(x, y, smth_spn,'sgolay');

    figure
    scatter(x, y, 'filled');
    hold on;
    scatter(x, yy, 'filled');
    hold off
 title(['Spec' num2str(DataNum)])
    %% find peaks
    min_dis = MeasurePara_Raman.minimum_peak_distance;

    % to check the result, execute only findpeaks(yy, x, 'MinPeakDistance', min_dis)
    [pks,locs, w, prom] = findpeaks(yy, x, 'MinPeakDistance', min_dis,'NPeaks',20 );
    [pks_neg,locs_neg, w_neg, prom_neg] = findpeaks(-yy, x, 'MinPeakDistance', min_dis,'NPeaks' ,10);
    


    %% Take a derivative 
    min_dis_deri = MeasurePara_Raman.minimum_peak_distance_deri;
    min_dis_deri2 = MeasurePara_Raman.minimum_peak_distance_deri2;

    dydx = diff(yy(:))./diff(x(:));
    dydx(end+1) = dydx(end);

    figure(5)
    scatter(x, dydx, 'filled');
    hold on;

    d2ydx2 = diff(dydx(:)) ./ diff(x(:));
    d2ydx2(end+1) = dydx(end);
    
    scatter(x, d2ydx2, 'filled');
    hold off
    
    % if the graph is negative, the function could be changed to positive to
    % negative. check if it is necessary by just run only findpeaks(-dydx, x,
    % 'MinPeakDistance', min_dis_deri).
   
    [pks_dy,locs_dy, w_dy, prom_dy] = findpeaks(dydx, x, 'MinPeakDistance', min_dis_deri,'NPeaks',50 );
    [pks_dy_neg,locs_dy_neg, w_dy_neg, prom_dy_neg] = findpeaks(-dydx, x, 'MinPeakDistance', min_dis_deri,'NPeaks' ,50);
    [pks_d2y,locs_d2y, w_d2y, prom_d2y] = findpeaks(d2ydx2, x, 'MinPeakDistance', min_dis_deri2,'NPeaks' ,50);
    [pks_d2y_neg,locs_d2y_neg, w_d2y_neg, prom_d2y_neg] = findpeaks(-d2ydx2, x, 'MinPeakDistance', min_dis_deri2,'NPeaks' ,50);
   
    figure(6)
    scatter(x, yy, 'filled','MarkerEdgeColor',[0 0.4470 0.7410])
    hold on
    scatter(locs,1.03*pks,100,'v','LineWidth',1.5,'MarkerEdgeColor',[0.8500 0.3250 0.0980])
    
     y_at_pks_neg = zeros(length(pks_neg),1);
    for j = 1:length(pks_neg)
        tmp_idx = find(x == locs_neg(j));
        y_at_pks_neg(j) = 0.97 * yy(tmp_idx);
    end
    scatter(locs_neg,y_at_pks_neg,100,'^','LineWidth',1.5,'MarkerEdgeColor',[0.8500 0.3250 0.0980])
    
    y_at_pks_dy = zeros(length(pks_dy),1);
    for j = 1:length(pks_dy)
        tmp_idx = find(x == locs_dy(j));
        y_at_pks_dy(j) = 1.03 * yy(tmp_idx);
    end
    scatter(locs_dy,y_at_pks_dy,100,'v','LineWidth',1.5,'MarkerEdgeColor',[0.9290 0.6940 0.1250])
    
        y_at_pks_dy_neg = zeros(length(pks_dy_neg),1);
    for j = 1:length(pks_dy_neg)
        tmp_idx = find(x == locs_dy_neg(j));
        y_at_pks_dy_neg(j) = 0.97 * yy(tmp_idx);
    end
    scatter(locs_dy_neg,y_at_pks_dy_neg,100,'^','LineWidth',1.5,'MarkerEdgeColor',[0.9290 0.6940 0.1250])
    
    y_at_pks_d2y = zeros(length(pks_d2y),1);
    for j = 1:length(pks_d2y)
        tmp_idx = find(x == locs_d2y(j));
        y_at_pks_d2y(j) = 1.03 * yy(tmp_idx);
    end
    scatter(locs_d2y,y_at_pks_d2y,100,'v','LineWidth',1.5,'MarkerEdgeColor',[0.4940 0.1840 0.5560])
    
        y_at_pks_d2y_neg = zeros(length(pks_d2y_neg),1);
    for j = 1:length(pks_d2y_neg)
        tmp_idx = find(x == locs_d2y_neg(j));
        y_at_pks_d2y_neg(j) = 0.97 * yy(tmp_idx);
    end
    scatter(locs_d2y_neg,y_at_pks_d2y_neg,100,'^','LineWidth',1.5,'MarkerEdgeColor',[0.4940 0.1840 0.5560])
    
    hold off
    legend('spec','pks','pks neg','pks dy','pks dy neg','pks d2y','pks d2y neg')
    title(['Spec' num2str(DataNum)])
%% Collecting feature values
features = [];
for i=1:length(list)
    i
    tst_dat(:,1) = MeasureData.x; tst_dat(:,2) = MeasureData.Data(:,i);
    
   
    x = tst_dat(tot_idx,1); y = tst_dat(tot_idx,2); y_raw = y;
    
     % background removal  -------comment out if it is UV-Vis.--------
     if MeasureData.Name ~= "UV/Vis"
         
         [BackCorrSpec] = msbackadj(x,y,'WindowSize',WindowSize,'StepSize',StepSize,.....
             'RegressionMethod','spline','EstimationMethod', 'quantile', 'QuantileValue', QuantileValue,'SmoothMethod', 'loess');
         nan_idx = logical(sum(isnan(BackCorrSpec),2));
         BackCorrSpec(nan_idx) = [];
         y = BackCorrSpec;
         x(nan_idx) = [];
     end

     % smoothing
     yy = smooth(x, y, smth_spn);
    
    % find peaks
    [pks,locs, w, prom] = findpeaks(yy, x, 'MinPeakDistance', min_dis,'NPeaks',2 );
    [pks_neg,locs_neg, w_neg, prom_neg] = findpeaks(-yy, x, 'MinPeakDistance', min_dis,'NPeaks' ,1);
    
    % first derivative
    dydx = diff(yy(:))./diff(x(:));
    dydx(end+1) = dydx(end);
    % second derivative
    d2ydx2 = diff(dydx(:))./diff(x(:));
    d2ydx2(end+1) = dydx(end);
    
    % derivative features
    [pks_dy,locs_dy, w_dy, prom_dy] = findpeaks(dydx, x, 'MinPeakDistance', min_dis_deri,'NPeaks',30 );
    [pks_dy_neg,locs_dy_neg, w_dy_neg, prom_dy_neg] = findpeaks(-dydx, x, 'MinPeakDistance', min_dis_deri,'NPeaks' ,30);
    [pks_d2y,locs_d2y, w_d2y, prom_d2y] = findpeaks(d2ydx2, x, 'MinPeakDistance', min_dis_deri2,'NPeaks' ,30);
    [pks_d2y_neg,locs_d2y_neg, w_d2y_neg, prom_d2y_neg] = findpeaks(-d2ydx2, x, 'MinPeakDistance', min_dis_deri2,'NPeaks' ,30);
    
    tmp = [pks', locs',w', prom', pks_neg',locs_neg', w_neg', prom_neg',....
        pks_dy', locs_dy', w_dy', prom_dy', pks_dy_neg', locs_dy_neg', w_dy_neg', prom_dy_neg',....
        pks_d2y', locs_d2y', w_d2y', prom_d2y', pks_d2y_neg', locs_d2y_neg', w_d2y_neg', prom_d2y_neg'];
    features = [features; tmp];
end

% %% Convert to a table
% FeatureData_Raman = array2table(features);
% FeatureData_Raman.Properties.VariableNames = {'pk1_int', 'pk2_int', 'pk1_loc', 'pk2_loc','pk1_w', 'pk2_w', 'pk1_prom', 'pk2_prom', ...
%     'vly1_int', 'vly1_loc', 'vly1_w', 'vly1_prom', ...
%     'pk1_dy_int', 'pk2_dy_int', 'pk3_dy_int','pk1_dy_loc', 'pk2_dy_loc','pk3_dy_loc','pk1_dy_w', 'pk2_dy_w', 'pk3_dy_w','pk1_dy_prom', 'pk2_dy_prom','pk3_dy_prom', ...
%     'vly1_dy_int', 'vly2_dy_int', 'vly3_dy_int','vly1_dy_loc', 'vly2_dy_loc','vly3_dy_loc','vly1_dy_w', 'vly2_dy_w', 'vly3_dy_w','vly1_dy_prom', 'vly2_dy_prom','vly3_dy_prom', ...
%     'pk1_d2y_int', 'pk2_d2y_int', 'pk3_d2y_int','pk1_d2y_loc', 'pk2_d2y_loc','pk3_d2y_loc','pk1_d2y_w', 'pk2_d2y_w', 'pk3_d2y_w','pk1_d2y_prom', 'pk2_d2y_prom','pk3_d2y_prom', ...
%     'vly1_d2y_int', 'vly2_d2y_int', 'vly3_d2y_int','vly1_d2y_loc', 'vly2_d2y_loc','vly3_d2y_loc','vly1_d2y_w', 'vly2_d2y_w', 'vly3_d2y_w','vly1_d2y_prom', 'vly2_d2y_prom','vly3_d2y_prom'};
% Target = randn(25,1); % input the target value here.
% tmpTbl = array2table(Target, 'VariableNames', "PEC");
% FeatureData = [FeatureData, tmpTbl];







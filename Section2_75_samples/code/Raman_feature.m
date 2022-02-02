%% Reading .CSV dataclea
list = dir('*.csv');

% addpath 'C:\Users\ideta\OneDrive - 中央大学';
load('MeasurePara');

%% Raman Opts
% インポート オプションの設定およびデータのインポート
opts = delimitedTextImportOptions("NumVariables", 3);

% 範囲と区切り記号の指定
opts.DataLines = [10, Inf];
opts.Delimiter = ",";

% 列名と型の指定
opts.VariableNames = ["DataCount2048", "Var2", "VarName3"];
opts.SelectedVariableNames = ["DataCount2048", "VarName3"];
opts.VariableTypes = ["double", "string", "double"];

% ファイル レベルのプロパティを指定
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% 変数プロパティを指定
opts = setvaropts(opts, "Var2", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var2", "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["DataCount2048", "VarName3"], "ThousandsSeparator", ",");
figure

for i = 1 : length(list)
    
    data(:,:,i) = readmatrix(list(i).name, opts); % read raw data
    plot(data(:,1,i),data(:,2,i))
    hold on
end
hold off
%%

MeasureData.x = data(:, 1, 1);
for i = 1: size(data, 3)
    MeasureData.Data(:, i) = data(:, 2, i);
end
MeasureData.Name = "Raman_spectrum";
MeasureData.xLabel = "Raman_Shift";

%% Parameters and conditions check using a single data

   %%
   % tst_dat = readmatrix(list(1).name, opts);
   tst_dat(:, 1) = MeasureData.x; tst_dat(:, 2) = MeasureData.Data(:, 56); % check by changing
   figure(1)
   % scatter(tst_dat(:, 1), tst_dat(:, 2), 'filled');←other representation
   plot(tst_dat(:, 1), tst_dat(:, 2));
   
   %% Range selection
   RamanShift_region = MeasurePara.Raman.RamanShift_region;
   shift_idx = find(tst_dat(:, 1) >= RamanShift_region(1) & tst_dat(:, 1) <= RamanShift_region(2));% select the raman shift region that I wanna use.
   
   figure(2)                                     
   % scatter(tst_dat(tot_idx,1), tst_dat(tot_idx,2), 'filled');
   plot(tst_dat(shift_idx, 1), tst_dat(shift_idx, 2));
   x = tst_dat(shift_idx, 1); y = tst_dat(shift_idx, 2);
   
  %% background elmination
   
  WindowSize = MeasurePara.Raman.WindowSize;
  StepSize = MeasurePara.Raman.StepSize;
  QuantileValue = MeasurePara.Raman.QuantileValue;
  
  [BackCorrSpec] = msbackadj(x, y, 'WindowSize', WindowSize, 'StepSize', StepSize,...
      'RegressionMethod', 'spline', 'EstimationMethod', 'quantile', 'QuantileValue', QuantileValue, 'SmoothMethod', 'loess');
 

  % To avoid error,remove Nan_data(not a number)
  Nan_idx = logical(sum(isnan(BackCorrSpec),2));
  BackCorrSpec(Nan_idx) = [ ];
  x(Nan_idx) = [ ];

  figure(3)
  scatter(x, y, 'filled');
  hold on
  scatter(x, BackCorrSpec, 'filled')
  hold off
  y = BackCorrSpec;
  
  %% take the smoothing within divided region
  xx = [ ];
  yy = [ ];
  dydx = [ ];
  dydx_smooth = [ ];
  d2ydx2 = [ ];
  peakReg = MeasurePara.Raman.Raman_Shift_Region;
  smoothing = MeasurePara.Raman.smoothing_span;
 for i = 1 : length(peakReg)
     % smoothing operation within each region
     peakRestr = peakReg(i, :);
     waveRestr = find((x >= peakRestr(1)) & (x <= peakRestr(2)));
     Cut_x = x(waveRestr);
     Cut_y = y(waveRestr);
     smoothing_span = smoothing(i);
     Cut_yy = smooth(Cut_x, Cut_y, smoothing_span, "sgolay");
%      % First derivative the data within each region
%      dydx_Cut = (diff(Cut_yy(:))./diff(Cut_x(:)));
%      dydx_Cut(end+1) = dydx_Cut(end);% fit the matrix size
%      dydx_smooth_Cut = smooth(Cut_x, dydx_Cut, smoothing_span, "sgolay");
%      % Second derivative the data within each region
%      d2ydx2_Cut = (diff(dydx_Cut(:))./diff(Cut_x(:)));
%      d2ydx2_Cut(end+1) = d2ydx2_Cut(end);
     % The data of the whole area is combined into one
     xx = [xx; Cut_x];
     yy = [yy; Cut_yy];
%      dydx = [dydx; dydx_Cut];
%      dydx_smooth = [dydx_smooth; dydx_smooth_Cut];
     
 end
 % plot the each data for confirming shape of figure
      figure(4)
     plot(x, y);
     hold on
     plot(xx, yy);
     
%      figure(5)
%      plot(xx, dydx);
%      hold on 
%      plot(xx, dydx_smooth);
%      plot(xx, d2ydx2);
%      hold off
 
%% Finding peaks operation
pks_posi = [ ];
locs_posi = [ ];
pks_neg = [ ];
locs_neg = [ ];
min_Dis = MeasurePara.Raman.minimum_peak_distance;
min_Dis_neg = MeasurePara.Raman.minimum_negativepeak_distance;
find_pks_reg = MeasurePara.Raman.find_peaks_region;
find_neg_pks_reg = MeasurePara.Raman.find_negativepeaks_region;
 % Positive peaks
 for j = 1 : length(find_pks_reg)
     minDis = min_Dis(j);
     peakfind = find_pks_reg(j, :);
     posi_waverange = find((xx >= peakfind(1)) & (xx <= peakfind(2)));
     Posi_x = xx(posi_waverange);
     Posi_yy = yy(posi_waverange);
       [pks_posi_Cut, locs_posi_Cut, w_posi_Cut, ~] = findpeaks(Posi_yy, Posi_x, 'MinPeakDistance', minDis);
     pks_posi= [pks_posi; pks_posi_Cut];
     locs_posi = [locs_posi; locs_posi_Cut];
 end
     figure(4)
     scatter(locs_posi, pks_posi, 'filled','^');
 % Negative peaks
 for k = 1 : length(find_neg_pks_reg)
     minDis_neg = min_Dis_neg(k);
     neg_peakfind = find_neg_pks_reg(k, :);
     neg_waverange = find((xx >= neg_peakfind(1)) & (xx <= neg_peakfind(2)));
     neg_x = xx(neg_waverange);
     neg_yy = yy(neg_waverange);
       [pks_neg_Cut, locs_neg_Cut, w_neg_Cut, ~] = findpeaks(-neg_yy, neg_x, 'MinPeakDistance', minDis_neg);
     pks_neg = [pks_neg; pks_neg_Cut];
     locs_neg = [locs_neg; locs_neg_Cut];
 end
     figure(4)
     scatter(locs_neg, -pks_neg, 'filled','v');
 
 %% This operation find sholder peaks
   min_dis_deri = MeasurePara.Raman.minimum_peak_distance_derivative;
   shol_Reg = MeasurePara.Raman.sholder_region;
   wave_shol = find((xx >= shol_Reg(1)) & (xx <= shol_Reg(2)));
   sholder_smoothing_span = MeasurePara.Raman.Sholder_smoothing_span;
   shol_x = xx(wave_shol);
   shol_yy = yy(wave_shol);
   shol_dydx = (diff(shol_yy(:))./diff(shol_x(:)));
   shol_dydx(end+1) = shol_dydx(end);% fit the matrix size
   shol_dydx_smooth = smooth(shol_x, shol_dydx, sholder_smoothing_span);
   % Find the derivative peaks
   [pks_dy_ori, locs_dy, w_dy, ~] = findpeaks(shol_dydx_smooth, shol_x, 'MinPeakDistance', min_dis_deri);
   locs_idx = find(xx == locs_dy);
   Pks_dy = yy(locs_idx);
   
   
   figure(4)
   scatter(locs_dy, Pks_dy, 'filled', 'd');
   
   figure(5)
   plot(shol_x,shol_dydx,'O-')
   hold on
   plot(shol_x,shol_dydx_smooth,'O-')
   scatter(xx(locs_idx),pks_dy_ori)
hold off
 %% Collecting feature values
 features = [ ];figure(6);hold on
 for i = 1 : length(list)
     tst_dat(:, 1) = MeasureData.x; tst_dat(:, 2) = MeasureData.Data(:, i);
     x = tst_dat(shift_idx,1); y = tst_dat(shift_idx,2); y_raw = y;
     i
     % background removal 
     [BackCorrSpec] = msbackadj(x, y, 'WindowSize', WindowSize, 'StepSize', StepSize,...
      'RegressionMethod', 'spline', 'EstimationMethod', 'quantile', 'QuantileValue', QuantileValue, 'SmoothMethod', 'loess');
  
     Nan_idx = logical(sum(isnan(BackCorrSpec),2));
     BackCorrSpec(Nan_idx) = [ ];
     y = BackCorrSpec;
     x(Nan_idx) = [ ];
     
 % take the smoothing within divided region
    xx = [ ];
    yy = [ ];
    dydx = [ ];
    dydx_smooth = [ ];
    d2ydx2 = [ ];
    peakReg = MeasurePara.Raman.Raman_Shift_Region;
    smoothing = MeasurePara.Raman.smoothing_span;
   for  j = 1 : length(peakReg)
     % smoothing operation within each region
     peakRestr = peakReg(j, :);
     waveRestr = find((x >= peakRestr(1)) & (x <= peakRestr(2)));
     Cut_x = x(waveRestr);
     Cut_y = y(waveRestr);
     smoothing_span = smoothing(j);
     Cut_yy = smooth(Cut_x, Cut_y, smoothing_span, "sgolay");
     % First derivative the data within each region
     dydx_Cut = (diff(Cut_yy(:))./diff(Cut_x(:)));
     dydx_Cut(end+1) = dydx_Cut(end);% fit the matrix size
     dydx_smooth_Cut = smooth(Cut_x, dydx_Cut, smoothing_span, "sgolay");
     % Second derivative the data within each region
     d2ydx2_Cut = (diff(dydx_Cut(:))./diff(Cut_x(:)));
     d2ydx2_Cut(end+1) = d2ydx2_Cut(end);
     % The data of the whole area is combined into one
     xx = [xx; Cut_x];
     yy = [yy; Cut_yy];
     dydx = [dydx; dydx_Cut];
     dydx_smooth = [dydx_smooth; dydx_smooth_Cut];
     d2ydx2 = [d2ydx2; d2ydx2_Cut];  
   end
   plot(xx,yy)
  % Finding peaks operation
     pks_posi = [ ];
     locs_posi = [ ];
     pks_neg = [ ];
     locs_neg = [ ];
     min_Dis = MeasurePara.Raman.minimum_peak_distance;
     min_Dis_neg = MeasurePara.Raman.minimum_negativepeak_distance;
     find_pks_reg = MeasurePara.Raman.find_peaks_region;
     find_neg_pks_reg = MeasurePara.Raman.find_negativepeaks_region;
 % Positive peaks
    for k = 1 : length(find_pks_reg)
     minDis = min_Dis(k);
     peakfind = find_pks_reg(k, :);
     posi_waverange = find((xx >= peakfind(1)) & (xx <= peakfind(2)));
     Posi_x = xx(posi_waverange);
     Posi_yy = yy(posi_waverange);
       [pks_posi_Cut, locs_posi_Cut, w_posi_Cut, ~] = findpeaks(Posi_yy, Posi_x, 'MinPeakDistance', minDis);
     pks_posi= [pks_posi; pks_posi_Cut];
     locs_posi = [locs_posi; locs_posi_Cut];
    end
    
  % Negative peaks
    for l = 1 : length(find_neg_pks_reg)
     minDis_neg = min_Dis_neg(l);
     neg_peakfind = find_neg_pks_reg(l, :);
     neg_waverange = find((xx >= neg_peakfind(1)) & (xx <= neg_peakfind(2)));
     neg_x = xx(neg_waverange);
     neg_yy = yy(neg_waverange);
       [pks_neg_Cut, locs_neg_Cut, w_neg_Cut, ~] = findpeaks(-neg_yy, neg_x, 'MinPeakDistance', minDis_neg);
     pks_neg = [pks_neg; pks_neg_Cut];
     locs_neg = [locs_neg; locs_neg_Cut];
    end
    
   % This operation find sholder peaks
   min_dis_deri = MeasurePara.Raman.minimum_peak_distance_derivative;
   shol_Reg = MeasurePara.Raman.sholder_region;
   wave_shol = find((xx >= shol_Reg(1)) & (xx <= shol_Reg(2)));
   sholder_smoothing_span = MeasurePara.Raman.Sholder_smoothing_span;
   shol_x = xx(wave_shol);
   shol_yy = yy(wave_shol);
   shol_dydx = (diff(shol_yy(:))./diff(shol_x(:)));
   shol_dydx(end+1) = shol_dydx(end);% fit the matrix size
   shol_dydx_smooth = smooth(shol_x, shol_dydx, sholder_smoothing_span);
   % Find the derivative peaks
   [pks_dy, locs_dy, w_dy, ~] = findpeaks(shol_dydx_smooth, shol_x, 'MinPeakDistance', min_dis_deri);
   locs_idx = find(xx == locs_dy);
   Pks_dy = yy(locs_idx);
     
     tmp = [pks_posi(1 : 5)', locs_posi(1 : 5)', pks_neg(1 : 4)', locs_neg(1 : 4)', Pks_dy, locs_dy];
     features = [features; tmp];
 end
 hold off;adjfig;xlabel('Ramanshift [cm^-^1]');ylabel('Intensity [-]')
 %% Convert to a table
 FeatureData_Raman = array2table(features);
     
 FeatureData_Raman.Properties.VariableNames = {'pk1_posi_int', 'pk2_posi_int', 'pk3_posi_int', 'pk4_posi_int', 'pk5_posi_int',...
     'pk1_posi_loc', 'pk2_posi_loc', 'pk3_posi_loc', 'pk4_posi_loc', 'pk5_posi_loc', 'pk_neg1_int', 'pk_neg2_int', 'pk_neg3_int', 'pk_neg4_int', ...
     'pk_neg1_loc', 'pk_neg2_loc', 'pk_neg3_loc', 'pk_neg4_loc', 'pk_dy_int', 'pk_dy_loc'};

 for i =1:size(FeatureData_Raman,2)
    FeatureData_Raman.Properties.VariableNames{i} = ['Raman_', FeatureData_Raman.Properties.VariableNames{i}];
end
%% 
%  addpath 'C:\Users\ideta\OneDrive - 中央大学'
%  load 'Features'
 FeatureData_Raman(:,11:18) =[];
     
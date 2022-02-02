% Reading .txt data
list=dir('*.txt');

% addpath 'C:\Users\ideta\OneDrive - 中央大学';
load('MeasurePara.mat');

% option for Nagai data
opts1 = delimitedTextImportOptions("NumVariables", 2);

% 範囲と区切り記号の指定
opts1.DataLines = MeasurePara.UVvis.x_cut1;
opts1.Delimiter = "\t";

% 列名と型の指定
opts1.VariableNames = ["Date", "VarName2"];
opts1.VariableTypes = ["double", "double"];

% ファイル レベルのプロパティを指定
opts1.ExtraColumnsRule = "ignore";
opts1.EmptyLineRule = "read";

% option for Idei & Chen data
opts2 = delimitedTextImportOptions("NumVariables", 2);

% 範囲と区切り記号の指定
opts2.DataLines = MeasurePara.UVvis.x_cut2;
opts2.Delimiter = "\t";

% 列名と型の指定
opts2.VariableNames = ["Date", "VarName2"];
opts2.VariableTypes = ["double", "double"];

% ファイル レベルのプロパティを指定
opts2.ExtraColumnsRule = "ignore";
opts2.EmptyLineRule = "read";
figure
for i=1:length(list)
    i
    if i <= 28  % Nagai samples
        opts = opts1;
    else  % Idei & Chen samples
        opts = opts2;
    end
    
    data(:,:,i) = readmatrix(list(i).name, opts); % raw data
    plot(data(:,1,i),data(:,2,i))
    hold on
end

hold off
%%


for i = 1 : size(data,3)
    MeasureData.x(:,i) = data(:, 1, i);
    MeasureData.Data(:, i) = data(:, 2, i);
end
MeasureData.Name = "UV/Vis";
MeasureData.xLabel = "Wavelength";

%% Parameters and conditions check using a single data
 %%
    % tst_dat = readmatrix(list(1).name, opts); 
    DataNum =17;
    tst_dat(:, 1) = MeasureData.x(:, DataNum); tst_dat(:, 2) = MeasureData.Data(:, DataNum); % to check by changing
    figure(1)
    % scatter(tst_dat(:,1), tst_dat(:,2), 'filled');
    scatter(tst_dat(:,1), tst_dat(:,2),'filled');
%% Range selection
    Wavelength_region = MeasurePara.UVvis.Wavelength_region;
    Removal_region1 = MeasurePara.UVvis.Removal_wavelength_region1;
    Removal_region2 = MeasurePara.UVvis.Removal_wavelength_region2;
    wave_idx = find((tst_dat(:,1) >= Wavelength_region(1)) & (tst_dat(:,1) <= Wavelength_region(2)));
    rem1_idx = find((tst_dat(:,1) < Removal_region1(1)) | (tst_dat(:,1) > Removal_region1(2))); 
    [tmp_idx, ~] = intersect(wave_idx, rem1_idx); 
    rem2_idx = find((tst_dat(:,1) < Removal_region2(1)) | (tst_dat(:,1) > Removal_region2(2)));
    [ovlp_idx, ~] = intersect(tmp_idx, rem2_idx); 
    
    figure(2)
    % scatter(tst_dat(ovlp_idx,1), tst_dat(ovlp_idx,2), 'filled');
    scatter(tst_dat(ovlp_idx,1), tst_dat(ovlp_idx,2),'filled');
    x = tst_dat(ovlp_idx,1); y = tst_dat(ovlp_idx,2);
    
    %% 
    % remove Nan(not a number)to avoid error
    nan_idx = logical(sum(isnan(y),2));
    y(nan_idx) = [];
    x(nan_idx) = [];
    
    figure(3)
    scatter(x,y,'filled')
    
%% Divide the wavelength region and smoothing/finding peak within each region
   %% Divide the wavelength region and Smoothing curve
   % In this case, I divided the wavelength region into 3 parts from 350nm
   % to 480nm, 480nm to 640nm and 640nm to 800nm
   Div_reg = MeasurePara.UVvis.division_region;
   smoothing_span = MeasurePara.UVvis.smoothing_span;
   smooth_x = [ ];
   smooth_y = [ ];
   smooth_yy = [ ];
   
   for i = 1 : length(Div_reg)
       div_reg = Div_reg(i, :);
       smooth_span = smoothing_span(i, :);
       Wave_reg = find((x >= div_reg(1)) & (x <= div_reg(2)));
       Cut_x = x(Wave_reg);
       Cut_y = y(Wave_reg);
       Cut_yy = smooth(Cut_x, Cut_y, smooth_span);
       smooth_x = [smooth_x; Cut_x];
       smooth_y = [smooth_y; Cut_y];
       smooth_yy = [smooth_yy; Cut_yy];
   end
   x = smooth_x; y = smooth_y; yy = smooth_yy;
    figure(4)
    scatter(x, y, 'filled');
%     plot(x, y);
    hold on;
    scatter(x, yy, 'filled');
%     plot(x, yy);
   
%% find maximum absorption point within certain region
  % decide the region which contain hematite peaks
  peak_reg = MeasurePara.UVvis.find_peaks_region;
  fdpk_reg = find((x >= peak_reg(1)) & (x <= peak_reg(2)));
  x_Cut = x(fdpk_reg);
  yy_Cut = yy(fdpk_reg);
  % find the maximum value and specify that index 
  yy_max = max(yy_Cut);
  max_ind = find(yy_Cut == yy_max);
  x_max = x_Cut(max_ind);
  
  pks = yy_max; locs = x_max;
  figure(4)
  scatter(locs, pks, '^', 'filled');
    
%% Derivative and find Sholder information
  % take a first derivative and smoothing derivation shapes for finding sholder peaks
  % Cut the wavelength region to derivate data
  shol_reg = MeasurePara.UVvis.sholder_peaks_region;
  
  
  shol_region = find((x >= shol_reg(1)) & (x <= shol_reg(2)));
  shol_x = x(shol_region);
  shol_yy = yy(shol_region);
  % first derivative and smoothing
  shol_smooth = MeasurePara.UVvis.sholder_smoothing_span;
  
  dydx = (diff(shol_yy(:))./diff(shol_x(:)));
  dydx(end+1) = dydx(end);
  smth_dydx = smooth(shol_x, dydx, shol_smooth);
   d2ydx2 = diff(smth_dydx(:)) ./ diff(shol_x(:));
    d2ydx2(end+1) = dydx(end);
  figure(5)
 plot(shol_x, dydx,'-o');
  hold on 
  plot(shol_x, smth_dydx,'-o');
plot(shol_x,d2ydx2)
  min_pks_dis_posi = MeasurePara.UVvis.minmum_peak_distance_posi; 
  min_pks_dis_neg = MeasurePara.UVvis.minmum_peak_distance_neg;
 
   
  % finding peaks operation within each region
  % positive derivative peaks
 [pks_posi_dy, locs_posi_dy, w_posi_dy, ~] = findpeaks(smth_dydx, shol_x, 'MinPeakDistance',  min_pks_dis_posi);

  % negative derivative peaks
 
  [pks_neg_dy, locs_neg_dy, w_neg_dy, ~] = findpeaks(-smth_dydx, shol_x, 'MinPeakDistance',  min_pks_dis_neg);
   
  figure(5)
  scatter(locs_posi_dy, pks_posi_dy, 'filled', 'b');
  scatter(locs_neg_dy, -pks_neg_dy, 'filled');

  % plot the peaks mark into figure(4)
  posipks_ind = find(x == locs_posi_dy);
  negpks_ind = find(x == locs_neg_dy);
  pks_posi_dy = yy(posipks_ind);
  pks_neg_dy = yy(negpks_ind);
    
  figure(4)
  scatter(locs_posi_dy, pks_posi_dy, 'filled', 'b');
  scatter(locs_neg_dy, pks_neg_dy, 'filled','y');
  
  figure(5)
  hold off
  %% Take average form 640nm to 800nm
    ave_region = MeasurePara.UVvis.average_region;
    ave_ind = find((x >= ave_region(1)) & (x <= ave_region(2)));
    ave_reg_x = x(ave_ind);
    ave_reg_abs = yy(ave_ind);
    average_locs = mean(ave_reg_x);
    average_abs = mean(ave_reg_abs);
    
    figure(4)
    scatter(average_locs, average_abs, 'filled', 's', 'r');
    hold off
 %% Collecting feature values       
 features = [ ];
 for i = 1 : length(list)
     i
   % Range selection
    tst_dat(:, 1) = MeasureData.x(:, i); tst_dat(:, 2) = MeasureData.Data(:, i);
    Wavelength_region = MeasurePara.UVvis.Wavelength_region;
    Removal_region1 = MeasurePara.UVvis.Removal_wavelength_region1;
    Removal_region2 = MeasurePara.UVvis.Removal_wavelength_region2;
    wave_idx = find((tst_dat(:,1) >= Wavelength_region(1)) & (tst_dat(:,1) <= Wavelength_region(2)));
    rem1_idx = find((tst_dat(:,1) < Removal_region1(1)) | (tst_dat(:,1) > Removal_region1(2))); 
    [tmp_idx, ~] = intersect(wave_idx, rem1_idx); 
    rem2_idx = find((tst_dat(:,1) < Removal_region2(1)) | (tst_dat(:,1) > Removal_region2(2)));
    [ovlp_idx, ~] = intersect(tmp_idx, rem2_idx); 
    
    x = tst_dat(ovlp_idx,1); y = tst_dat(ovlp_idx,2);
    
  % 
    % remove Nan(not a number)to avoid error
    nan_idx = logical(sum(isnan(y),2));
    y(nan_idx) = [];
    x(nan_idx) = [];
    
 % Divide the wavelength region and smoothing/finding peak within each region
  % Divide the wavelength region and Smoothing curve
   % In this case, I divided the wavelength region into 3 parts from 350nm
   % to 480nm, 480nm to 640nm and 640nm to 800nm
   Div_reg = MeasurePara.UVvis.division_region;
   smoothing_span = MeasurePara.UVvis.smoothing_span;
   smooth_x = [ ];
   smooth_y = [ ];
   smooth_yy = [ ];
   
   for j = 1 : length(Div_reg)
       div_reg = Div_reg(j, :);
       smooth_span = smoothing_span(j, :);
       Wave_reg = find((x >= div_reg(1)) & (x <= div_reg(2)));
       Cut_x = x(Wave_reg);
       Cut_y = y(Wave_reg);
       Cut_yy = smooth(Cut_x, Cut_y, smooth_span, 'sgolay');
       smooth_x = [smooth_x; Cut_x];
       smooth_y = [smooth_y; Cut_y];
       smooth_yy = [smooth_yy; Cut_yy];
   end
   x = smooth_x; y = smooth_y; yy = smooth_yy;
   
 % find maximum absorption point within certain region
  % decide the region which contain hematite peaks
  peak_reg = MeasurePara.UVvis.find_peaks_region;
  fdpk_reg = find((x >= peak_reg(1)) & (x <= peak_reg(2)));
  x_Cut = x(fdpk_reg);
  yy_Cut = yy(fdpk_reg);
  % find the maximum value and specify that index 
  yy_max = max(yy_Cut);
  max_ind = find(yy_Cut == yy_max);
  x_max = x_Cut(max_ind);
  
  pks = yy_max; locs = x_max;
    
 % Derivative and find Sholder information
  % take a first derivative and smoothing derivation shapes for finding sholder peaks
  % Cut the wavelength region to derivate data
  shol_reg = MeasurePara.UVvis.sholder_peaks_region;
  shol_region = find((x >= shol_reg(1)) & (x <= shol_reg(2)));
  shol_x = x(shol_region);
  shol_yy = yy(shol_region);
  % first derivative and smoothing
  shol_smooth = MeasurePara.UVvis.sholder_smoothing_span;
  
  dydx = (diff(shol_yy(:))./diff(shol_x(:)));
  dydx(end+1) = dydx(end);
  smth_dydx = smooth(shol_x, dydx, shol_smooth);
  
  min_pks_dis_posi = MeasurePara.UVvis.minmum_peak_distance_posi; 
  min_pks_dis_neg = MeasurePara.UVvis.minmum_peak_distance_neg;
 
  % finding peaks operation within each region
  % positive derivative peaks
  [pks_posi_dy, locs_posi_dy, w_posi_dy, ~] = findpeaks(smth_dydx, shol_x, 'MinPeakDistance',  min_pks_dis_posi);
  % negative derivative peaks
  [pks_neg_dy, locs_neg_dy, w_neg_dy, ~] = findpeaks(-smth_dydx, shol_x, 'MinPeakDistance',  min_pks_dis_neg);

  % plot the peaks mark into figure(4)
  posipks_ind = find(x == locs_posi_dy);
  negpks_ind = find(x == locs_neg_dy);
  pks_posi_dy = yy(posipks_ind);
  pks_neg_dy = yy(negpks_ind);
  
  % Take average form 640nm to 800nm
    ave_region = MeasurePara.UVvis.average_region;
    ave_ind = find((x >= ave_region(1)) & (x <= ave_region(2)));
    ave_reg_x = x(ave_ind);
    ave_reg_abs = yy(ave_ind);
    average_locs = mean(ave_reg_x);
    average_abs = mean(ave_reg_abs);
    
 tmp = [pks, locs, pks_posi_dy, locs_posi_dy, pks_neg_dy, locs_neg_dy, average_abs];
 features = [features; tmp];
 end
 
 %% Convert to a table
FeatureData_UV_Vis = array2table(features);
FeatureData_UV_Vis.Properties.VariableNames = {'pks_abs', 'pks_locs', 'shol_posi_abs', 'shol_posi_locs',...
    'shol_neg_abs', 'shol_neg_locs', 'average_abs'};

for i =1:size(FeatureData_UV_Vis,2)
    FeatureData_UV_Vis.Properties.VariableNames{i} = ['UV_Vis_', FeatureData_UV_Vis.Properties.VariableNames{i}];
end
%% 
%  addpath 'C:\Users\ideta\OneDrive - 中央大学'
%  load 'Features'
%  Features.UVvis.No1_103 = FeatureData_UVVis;
%     
    
    
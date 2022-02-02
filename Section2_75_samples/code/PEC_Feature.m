


list1 = dir('*txt');

% インポート オプションの設定およびデータのインポート
opts1 = delimitedTextImportOptions("NumVariables", 2);

% 範囲と区切り記号の指定
DataLines = [20, Inf];

opts1.Delimiter = " ";

% 列名と型の指定
opts1.VariableNames = ["PotentialV", "CurrentA"];
opts1.VariableTypes = ["double", "double"];

% ファイル レベルのプロパティを指定
opts1.ExtraColumnsRule = "ignore";
opts1.EmptyLineRule = "read";
opts1.ConsecutiveDelimitersRule = "join";
opts1.LeadingDelimitersRule = "ignore";

% 変数プロパティを指定
opts1 = setvaropts(opts1, "PotentialV", "TrimNonNumeric", true);
opts1 = setvaropts(opts1, "PotentialV", "ThousandsSeparator", ",");


% データのインポート
for i = 1:length(list1)
    i

   
    data_tmp = readtable(list1(i).name, opts1);
    data1_tmp = table2array(data_tmp);
    nan_idx = (isnan(data1_tmp));
    data1_tmp(nan_idx(:,1),:) =[];
    data1(:,:,i) = data1_tmp;
    clear data_tmp 
end

vol1 = data1(:,1,:) + 0.059 * MeasurePara_PEC.pH + 0.1976;
vol1 = squeeze(vol1);
pec1 = 10^3 * data1(:,2,:);
pec1 = squeeze(pec1)./(MeasurePara_PEC.surface_area)';

for i =1:length(list1)
    [~,idx_performance] = min(abs(1.23-vol1(:,i)));
    Target1(i,1) = pec1(idx_performance,i);
end
    
list2 = dir('*csv');

% インポート オプションの設定およびデータのインポート
opts2 = delimitedTextImportOptions("NumVariables", 5);

% 範囲と区切り記号の指定
opts2.DataLines = [5, Inf];
opts2.Delimiter = ",";

% 列名と型の指定
opts2.VariableNames = ["Var1", "Var2", "OutputVoltageV", "CurrentA", "Var5"];
opts2.SelectedVariableNames = ["OutputVoltageV", "CurrentA"];
opts2.VariableTypes = ["string", "string", "double", "double", "string"];

% ファイル レベルのプロパティを指定
opts2.ExtraColumnsRule = "ignore";
opts2.EmptyLineRule = "read";

% 変数プロパティを指定
opts2 = setvaropts(opts2, ["Var1", "Var2", "Var5"], "WhitespaceRule", "preserve");
opts2 = setvaropts(opts2, ["Var1", "Var2", "Var5"], "EmptyFieldRule", "auto");


for i = length(list1) + 1:length(list1) +length(list2)
    i
    data_tmp = readtable(list2(i-28).name, opts2);
    
    % 出力型への変換
    data2_tmp = table2array(data_tmp);
data2_tmp = data2_tmp(end-672:end,:);
    data2(:,:,i-28) = data2_tmp;
    clear data_tmp 
end

vol2 = data2(:,1,:) + 0.059 * MeasurePara_PEC.pH + 0.1976;
vol2 = squeeze(vol2);
pec2 = 10^3 * data2(:,2,:);
pec2 = squeeze(pec2);

for i = length(list1) + 1:length(list1) +length(list2)
    
    [~,idx_performance] = min(abs(1.23-vol2(:,i-28)));
    Target2(i-28,1) = pec2(idx_performance,i-28);
end
    

Target = array2table(vertcat(Target1,Target2));
 Target.Properties.VariableNames{1} = 'PEC';


%% import UV/Vis
cd('UV_Vis')

[x,y,abs_files] = getuv_ref;

All_data_bare_hem.UV_vis.wavelength = x;
All_data_bare_hem.UV_vis.absorbance = y;
All_data_bare_hem.UV_vis.files = abs_files;

clear x y
cd('..//')

%% import Raman
cd('Raman\select')

[x,y,raman_files] = getraman_uec;

All_data_bare_hem.Raman.raman_shift = x;
All_data_bare_hem.Raman.intensity = y;
All_data_bare_hem.Raman.files = raman_files;

clear x y
cd('..//..//')

%% import PEIS
cd('PEIS')
[z_re,z_im,z_area,peis_files,imp_abs,freq,phase_shift] = getimp;

All_data_bare_hem.PEIS.freq = freq;
All_data_bare_hem.PEIS.z_re = z_re;
All_data_bare_hem.PEIS.z_im = z_im;
All_data_bare_hem.PEIS.imp = imp_abs;
All_data_bare_hem.PEIS.phase_shift = phase_shift;
All_data_bare_hem.PEIS.files = peis_files;

clearvars -except All_data_bare_hem  
cd('..//')
%% import XRD
cd("XRD")

[x,y,legend_txt,xrd_files] = getxrd;

All_data_bare_hem.XRD.diff_angle = x;
All_data_bare_hem.XRD.intensity = y;
All_data_bare_hem.XRD.files = xrd_files;

clear x y
cd('..//')


%% import target (PEC)

cd('PEC')

[x,y,iv_files] = getiv;

xlsm_files = dir('*.xlsm');

surface_area = zeros(length(xlsm_files),1);
for i =1:length(xlsm_files)

surface_area(i) = readmatrix(xlsm_files(i).name,'Range','F1:F1');

end
All_data_bare_hem.PEC.vol = x;
All_data_bare_hem.PEC.current = y;
All_data_bare_hem.PEC.surface_area = surface_area;
All_data_bare_hem.PEC.current_density = y./surface_area';
All_data_bare_hem.PEC.files = iv_files;

clearvars -except All_data_bare_hem  
cd('..//')




%% import functions 

% UV/Vis
function [x,y,abs_files] = getuv_ref(varargin)
%% control parameters
% Data files
ext='txt';     %データの拡張子
data_dir=[];   %データのディレクトリ
include_and=[];  % example include_and=[{'dark'},{'600'}] or [] ...both 600 and dark
exclude_and=[]; % example add=[{'dark'},{'600'}] or [] ...not both 600 and dark
include_or=[]; % example include_or=[{'dark'},{'600'}] or [] ..either 600 or dark
exclude_or=[];% example exclude_or=[{'dark'},{'600'}] or []..neither 600 nor dark
intervals=[]; % [] or [a,b].... Only a*x+b th data which satisfies above conditions are extracted  
move_data_dir=[];%'C:\Users\YuyaNagai\Desktop\Lab\Research\test'

% figures
font_name='arial';
label_y='Absorbance [-]';   %y軸タイトル
label_x='Wavelength [nm]';   %x軸タイトル
font_sz=20;   %図のフォントサイズ
line_width=1.5;   %プロット線の太さ
legend_box='off';   %凡例のボックスの表示('on')、非表示('off')
fig_copy='off';    %図のクリップボードへのコピーあり('on')、なし('off')

% save
save_dir=['C:\Users\YuyaNagai\Desktop\Lab\Research\test\test'];    %計算結果等を保存するディレクトリ
single_save_id='off';    %セーブあり('on')、なし('off')、上書きなし
multiple_save_id='off';    %連番をつけたうえでのセーブあり('on')、なし('off')
func_back_up='off';    %使用した関数のバックアップ(txt形式)あり('on')、なし('off')
multiple_func_back_up='off';    %連番での関数のバックアップ(txt形式)あり('on')、なし('off')

fig_file_name='raw_data'; %保存する図のタイトル
calc_file_name='raw_data'; %保存するデータのタイトル
save_variables={'x','y','abs_files'}; %保存する変数
file_count_method='%05i';  %連番の設定

% other

pH=13.61;

color_config=[        0         0    1.0000
    1.0000         0         0
    0    0.5000         0
    1.0000         0    1.0000
    0.5000         0    1.0000
    0    0.7500    1.0000
    1.0000    0.5000         0
    0.5000    0.7500         0
    1.0000    0.7500    0.7500
    0.5000    0.5000    1.0000
    0.7500    0.7500         0
    0.5000    0.5000    0.5000
    0.7500         0         0
    0         0         0];
% desing the color matrix. column corresponds to R, G, B





%% Setting
% data directly
oldfolder=cd;
if isempty(data_dir)~=1
    
    try
        cd(data_dir);
        
    catch
        mkdir(data_dir);
        cd(data_dir);
        
    end
else
    data_dir=cd;
    cd(data_dir)
end

func_dir=mfilename('fullpath');
func_name=mfilename;
previous_func_back_up=dir(['*',func_name,'*',]);
func_back_up_file_name=[func_name,'_back_up'];

num_ext=length(ext);

figure1 = figure('Color',[1 1 1],'visible','off');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
data_files_all=dir(['*',ext]);
data_files_all=data_files_all(~ismember({data_files_all.name},{previous_func_back_up.name}));
%% data search, filtering
% searching specific (and,or) data
if isempty(varargin)==1
    
    
    % and search
    % 1) include
    if isempty(include_and)==1
        idx_include_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(include_and)
            if i==1
                idx_include_and=ones(1,length(data_files_all));
            end
            include_and_tmp=['*',include_and{i},'*'];
            include_and_files=dir(include_and_tmp);
            idx_include_and_tmp=ismember({data_files_all.name},{include_and_files.name});
            idx_include_and=(idx_include_and.*idx_include_and_tmp);
            
        end
        
    end
    
    
    % 2) exclude
    if isempty(exclude_and)==1
        idx_exclude_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_and)
            if i==1
                idx_exclude_and=ones(1,length(data_files_all));
            end
            exclude_and_tmp=['*',exclude_and{i},'*'];
            exclude_and_files=dir(exclude_and_tmp);
            idx_exclude_and_tmp=ismember({data_files_all.name},{exclude_and_files.name});
            idx_exclude_and=(idx_exclude_and.*idx_exclude_and_tmp);
            
        end
        idx_exclude_and=~idx_exclude_and;
    end
    idx_and=((idx_include_and.*idx_exclude_and)>=1)';
    
    
    % or search
    % 1) include
    if isempty(include_or)==1
        idx_include_or=ones(1,length(data_files_all));
        
    else
        
        for i=1:length(include_or)
            if i==1
                idx_include_or=ones(1,length(data_files_all));
            end
            include_or_tmp=['*',include_or{i},'*'];
            include_or_files=dir(include_or_tmp);
            idx_include_or_tmp=~ismember({data_files_all.name},{include_or_files.name});
            idx_include_or=(idx_include_or.*idx_include_or_tmp);
            
        end
        idx_include_or=~idx_include_or;
    end
    
    % 2) exclude
    if isempty(exclude_or)==1
        idx_exclude_or=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_or)
            if i==1
                idx_exclude_or=logical(length(data_files_all));
            end
            exclude_or_tmp=['*',exclude_or{i},'*'];
            exclude_or_files=dir(exclude_or_tmp);
            idx_exclude_or_tmp=~ismember({data_files_all.name},{exclude_or_files.name});
            idx_exclude_or=(idx_exclude_or.*idx_exclude_or_tmp);
            
        end
    end
    
    idx_or=((idx_include_or.*idx_exclude_or)>=1)';
    
    
    idx2see=(idx_and.*idx_or)>=1;
    
    abs_files=data_files_all(idx2see);
    
    
else
    file_name=varargin;
    data_files_all(1:end)=[];
    abs_files=data_files_all;
    for j=1:length(file_name)
        
        abs_files(j)=dir(file_name{j});
    end
    
end

if isempty(intervals)~=1
    max_num_files=length(abs_files);
    
    counting_number=0;
    tmp_files=dir(ext);
    while max_num_files>=intervals(1)*counting_number+intervals(2)
        tmp_files(counting_number+1,1)=abs_files(intervals(1)*counting_number+intervals(2));
        counting_number=counting_number+1;
    end
    
    abs_files=tmp_files;
end
num_files=length(abs_files);

if num_files==0
    disp('No data satisfies the conditions you set')
    close gcf
    return
end
%%  データのインポート(以下個々のデータでコード生成して張り付け、ループにする)
for i=1:num_files
    %% 入力の取り扱い
    dataLines = [20, Inf];
    %     % dataLines が指定されていない場合、既定値を定義します
    %     if nargin < 2
    %         dataLines = [20, Inf];
    %     end
    %
  %% インポート オプションの設定およびデータのインポート
opts = delimitedTextImportOptions("NumVariables", 2);

% 範囲と区切り記号の指定
opts.DataLines = [16, Inf];
opts.Delimiter = "\t";

% 列名と型の指定
opts.VariableNames = ["Date", "VarName2"];
opts.VariableTypes = ["double", "double"];

% ファイル レベルのプロパティを指定
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% データのインポート
mydata{i} = readtable(abs_files(i).name, opts);

%% 出力型への変換
    raw_data{i} = table2array(mydata{i});
    
    x_cell_tmp=raw_data{i}(:,1);
    y_cell_tmp=raw_data{i}(:,2);
    x_cell_tmp(isnan(x_cell_tmp))=[];
    y_cell_tmp(isnan(y_cell_tmp))=[];
    x_cell{i}=x_cell_tmp;
    y_cell{i}=y_cell_tmp;
    data_length(i)=length(x_cell{i});
    %% グラフ描画
    plot(x_cell{i},y_cell{i},'LineWidth',1.5)
    
    % ファイル名から凡例を作成
    legend_tmp=abs_files(i).name;
    L{i}=find(legend_tmp=='_');
    legend_tmp(L{i})=' ';
    legend_tmp(end-num_ext:end)=[]; % ファイル名の拡張子部分の削除
    legend_txt{i}=legend_tmp;
    
 %% データ移動(コピー)
 
 if isempty(move_data_dir)==0
     
     try
         cd(move_data_dir)
         move_data_dir=cd;
         copyfile([abs_files(i).folder,'\',abs_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new','.',ext])
     catch
         cd(oldfolder)
         mkdir(move_data_dir);
         cd(move_data_dir);
         move_data_dir=cd;
         copyfile([abs_files(i).folder,'\',abs_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new.','.',ext])
     end
 end
end
%% グラフ及び出力値の調整
if length(unique(data_length))==1 && numel(data_length)>1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
    for l=1:num_files-1
        x_id(:,l)=x(:,l+1)-x(:,l);
    end
    if sum(abs(x_id(:)))==0
        x=x(:,1);
    end
    
elseif numel(data_length)==1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
else
    x=x_cell;
    y=y_cell;
    
end


% ylabel を作成
ylabel(label_y);
% xlabel を作成
xlabel(label_x);
box(axes1,'on');
% 残りの座標軸プロパティの設定
set(axes1,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
% legend を作成
if num_files<20
    legend1=legend(legend_txt,'Box',legend_box,'FontName',font_name);
end
hold off
colororder(color_config)
if strcmp(fig_copy,'on')==1  % クリップボードにコピーする
    % 現在の図の背景などを透明色に指定
set(gcf,'InvertHardcopy','off')
set(gcf,'Color','none');
set(gca,'Color','none');
% クリップボードにコピーして背景などを元に戻す
try
print('-clipboard','-dmeta')
catch
    print('-clipboard','-dbitmap')
end
set(gcf,'InvertHardcopy','on')
set(gcf,'Color',[1,1,1]);
set(gca,'Color','[1,1,1]');
end
set(gcf,'visible','on')
%%  データ及び図の保存
tr=strcmp(single_save_id,'on');
tr_multi=strcmp(multiple_save_id,'on');
tr_func_back_up=strcmp(func_back_up,'on');
tr_func_back_up_multi=strcmp(multiple_func_back_up,'on');
if tr==1 || tr_multi==1 || tr_func_back_up==1 || tr_func_back_up_multi==1
    if isempty(save_dir)==1
        save_dir=oldfolder;
        cd(save_dir)
    end
    try
        cd(save_dir);
        save_dir=cd;
    catch
        cd(oldfolder)
        mkdir(save_dir);
        cd(save_dir);
        save_dir=cd;
    end
    
    
    previous_fig_file=dir(['*',fig_file_name,'*.fig']);
    previous_calc_file=dir(['*',calc_file_name,'*.mat']);
    previous_func_back_up=dir(['*',func_name,'*',]);
    
    if tr_multi==1 || tr_func_back_up_multi==1
        
        
        file_count_fig=sprintf(file_count_method,length(previous_fig_file)+1);
        file_count_calc=sprintf(file_count_method,length(previous_calc_file)+1);
        file_count_func=sprintf(file_count_method,length(previous_func_back_up)+1);
        
        fig_file_name=[fig_file_name,'_',file_count_fig];
        calc_file_name=[calc_file_name,'_',file_count_calc];
        func_back_up_file_name=[func_back_up_file_name,'_',file_count_func];
    end
    
    
    
    
    calc_file_name=append("'",calc_file_name,"'");
    
    for k=1:length(save_variables)
        
        
        save_variables_tmp1=(save_variables{k});
        
        save_variables_tmp{k}=append(",","'",save_variables{k},"'");
        expression_assign=join(["assignin('caller'",save_variables_tmp{k},",",save_variables_tmp1,")"]);
        eval(expression_assign)
        
    end
    save_variables=join(string(save_variables_tmp));
    
    expression=join(["save(",calc_file_name,save_variables,")"]);
    % expression...save( 'raw_data' , 'x','y','iv_files' )
    
    if isempty(previous_fig_file)==1
        savefig(fig_file_name);
    end
    
    if isempty(previous_calc_file)==1
        evalin('caller',expression)  % save the variales by evalin
    end
    
    
end


if tr_func_back_up==1 || tr_func_back_up_multi==1
    copyfile([func_dir,'.m'],[save_dir,'\',func_back_up_file_name,'.txt'])
end



cd(oldfolder)





end

% XRD
function [x,y,legend_txt,xrd_files] = getxrd(varargin)
    %% control parameters
% Data files
ext='ras';     %データの拡張子
data_dir=[];   %データのディレクトリ
include_and=[];  % example include_and=[{'dark'},{'600'}] or [] ...both 600 and dark
exclude_and=[]; % example add=[{'dark'},{'600'}] or [] ...not both 600 and dark
include_or=[]; % example include_or=[{'dark'},{'600'}] or [] ..either 600 or dark
exclude_or=[];% example exclude_or=[{'dark'},{'600'}] or []..neither 600 nor dark
intervals=[]; % [] or [a,b].... Only a*x+b th data which satisfies above conditions are extracted  
move_data_dir=[];

% figures
font_name='arial';
label_y='Intensity [-]';   %y軸タイトル
label_x='2\theta [degree]';   %x軸タイトル
font_sz=20;   %図のフォントサイズ
line_width=1.5;   %プロット線の太さ
legend_box='off';   %凡例のボックスの表示('on')、非表示('off')
fig_copy='on';    %図のクリップボードへのコピーあり('on')、なし('off')

% save
save_dir=[];    %計算結果等を保存するディレクトリ
single_save_id='off';    %セーブあり('on')、なし('off')、上書きなし
multiple_save_id='off';    %連番をつけたうえでのセーブあり('on')、なし('off')
func_back_up='off';    %使用した関数のバックアップ(txt形式)あり('on')、なし('off')
multiple_func_back_up='off';    %連番での関数のバックアップ(txt形式)あり('on')、なし('off')

fig_file_name='raw_data'; %保存する図のタイトル
calc_file_name='raw_data'; %保存するデータのタイトル
save_variables={'x','y','xrd_files'}; %保存する変数
file_count_method='%05i';  %連番の設定

% other


color_config=[        0         0    1.0000
    1.0000         0         0
    0    0.5000         0
    1.0000         0    1.0000
    0.5000         0    1.0000
    0    0.7500    1.0000
    1.0000    0.5000         0
    0.5000    0.7500         0
    1.0000    0.7500    0.7500
    0.5000    0.5000    1.0000
    0.7500    0.7500         0
    0.5000    0.5000    0.5000
    0.7500         0         0
    0         0         0];
% desing the color matrix. column corresponds to R, G, B





%% Setting
% data directly
oldfolder=cd;
if isempty(data_dir)~=1
    
    try
        cd(data_dir);
        
    catch
        mkdir(data_dir);
        cd(data_dir);
        
    end
else
    data_dir=cd;
    cd(data_dir)
end

func_dir=mfilename('fullpath');
func_name=mfilename;
previous_func_back_up=dir(['*',func_name,'*',]);
func_back_up_file_name=[func_name,'_back_up'];

num_ext=length(ext);

figure1 = figure('Color',[1 1 1],'visible','off');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
data_files_all=dir(['*',ext]);
data_files_all=data_files_all(~ismember({data_files_all.name},{previous_func_back_up.name}));
%% data search, filtering
% searching specific (and,or) data
if isempty(varargin)==1
    
    
    % and search
    % 1) include
    if isempty(include_and)==1
        idx_include_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(include_and)
            if i==1
                idx_include_and=ones(1,length(data_files_all));
            end
            include_and_tmp=['*',include_and{i},'*'];
            include_and_files=dir(include_and_tmp);
            idx_include_and_tmp=ismember({data_files_all.name},{include_and_files.name});
            idx_include_and=(idx_include_and.*idx_include_and_tmp);
            
        end
        
    end
    
    
    % 2) exclude
    if isempty(exclude_and)==1
        idx_exclude_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_and)
            if i==1
                idx_exclude_and=ones(1,length(data_files_all));
            end
            exclude_and_tmp=['*',exclude_and{i},'*'];
            exclude_and_files=dir(exclude_and_tmp);
            idx_exclude_and_tmp=ismember({data_files_all.name},{exclude_and_files.name});
            idx_exclude_and=(idx_exclude_and.*idx_exclude_and_tmp);
            
        end
        idx_exclude_and=~idx_exclude_and;
    end
    idx_and=((idx_include_and.*idx_exclude_and)>=1)';
    
    
    % or search
    % 1) include
    if isempty(include_or)==1
        idx_include_or=ones(1,length(data_files_all));
        
    else
        
        for i=1:length(include_or)
            if i==1
                idx_include_or=ones(1,length(data_files_all));
            end
            include_or_tmp=['*',include_or{i},'*'];
            include_or_files=dir(include_or_tmp);
            idx_include_or_tmp=~ismember({data_files_all.name},{include_or_files.name});
            idx_include_or=(idx_include_or.*idx_include_or_tmp);
            
        end
        idx_include_or=~idx_include_or;
    end
    
    % 2) exclude
    if isempty(exclude_or)==1
        idx_exclude_or=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_or)
            if i==1
                idx_exclude_or=logical(length(data_files_all));
            end
            exclude_or_tmp=['*',exclude_or{i},'*'];
            exclude_or_files=dir(exclude_or_tmp);
            idx_exclude_or_tmp=~ismember({data_files_all.name},{exclude_or_files.name});
            idx_exclude_or=(idx_exclude_or.*idx_exclude_or_tmp);
            
        end
    end
    
    idx_or=((idx_include_or.*idx_exclude_or)>=1)';
    
    
    idx2see=(idx_and.*idx_or)>=1;
    
    xrd_files=data_files_all(idx2see);
    
    
else
    file_name=varargin;
    data_files_all(1:end)=[];
    xrd_files=data_files_all;
    for j=1:length(file_name)
        
        xrd_files(j)=dir(file_name{j});
    end
    
end

if isempty(intervals)~=1
    max_num_files=length(xrd_files);
    
    counting_number=0;
    tmp_files=dir(ext);
    while max_num_files>=intervals(1)*counting_number+intervals(2)
        tmp_files(counting_number+1,1)=xrd_files(intervals(1)*counting_number+intervals(2));
        counting_number=counting_number+1;
    end
    
    xrd_files=tmp_files;
end
num_files=length(xrd_files);

if num_files==0
    disp('No data satisfies the conditions you set')
    close gcf
    return
end
%%  データのインポート(以下個々のデータでコード生成して張り付け、ループにする)   
    for i=1:num_files
       
            



%IMPORTFILE2 テキスト ファイルからデータをインポート
%  UNTITLED = IMPORTFILE2(FILENAME) は既定の選択に関してテキスト ファイル FILENAME
%  からデータを読み取ります。  データを table として返します。
%
%  UNTITLED = IMPORTFILE2(FILE, DATALINES) はテキスト ファイル FILENAME
%  の指定された行区間のデータを読み取ります。DATALINES
%  を正の整数スカラーとして指定するか、行区間が不連続の場合は正の整数スカラーからなる N 行 2 列の配列として指定します。
%
%  例:
%  Untitled = importfile2("C:\Users\YuyaNagai\Desktop\Lab\Research\XRD\4mix\XRD\20200715_01_1.ras", [319, 3819]);
%
%  READTABLE も参照してください。
%
% MATLAB からの自動生成日: 2020/07/16 12:10:33

%% 入力の取り扱い

% dataLines が指定されていない場合、既定値を定義します

    dataLines = [319, Inf];


%% インポート オプションの設定およびデータのインポート
opts = delimitedTextImportOptions("NumVariables", 2);

% 範囲と区切り記号の指定
opts.DataLines = dataLines;
opts.Delimiter = " ";

% 列名と型の指定
opts.VariableNames = ["VarName1", "VarName2"];
opts.VariableTypes = ["string", "string"];

% ファイル レベルのプロパティを指定
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% 変数プロパティを指定
opts = setvaropts(opts, ["VarName1", "VarName2"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName1", "VarName2"], "EmptyFieldRule", "auto");

% データのインポート
tmp_data = readtable(xrd_files(i).name, opts);


   %% 出力型への変換
% if i>1
%     data{i}(:,1)=[];
% end

nodata_ind=ismissing(tmp_data(:,2));

tmp_data(nodata_ind,:)=[];

data{i} =str2double(table2array(tmp_data));

x_cell{i}=data{i}(:,1);
y_cell{i}=data{i}(:,2);
 data_length(i)=length(x_cell{i});

% all_data=cell2mat(data);
% x=all_data(:,1);
% y=all_data(:,2:end);

  %% グラフ描画
    plot(x_cell{i},y_cell{i},'LineWidth',1.5)
    
    % ファイル名から凡例を作成
    legend_tmp=xrd_files(i).name;
    L{i}=find(legend_tmp=='_');
    legend_tmp(L{i})=' ';
    legend_tmp(end-num_ext:end)=[]; % ファイル名の拡張子部分の削除
    legend_txt{i}=legend_tmp;

 %% データ移動(コピー)
 if isempty(move_data_dir)==0
     
     try
         cd(move_data_dir)
         move_data_dir=cd;
         copyfile([xrd_files(i).folder,'\',xrd_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new','.',ext])
     catch
         cd(oldfolder)
         mkdir(move_data_dir);
         cd(move_data_dir);
         move_data_dir=cd;
         copyfile([xrd_files(i).folder,'\',xrd_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new.','.',ext])
     end
 
 end
    end
 
    %% グラフ及び出力値の調整
if length(unique(data_length))==1 && numel(data_length)>1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
    for l=1:num_files-1
        x_id(:,l)=x(:,l+1)-x(:,l);
    end
    if sum(abs(x_id(:)))==0
        x=x(:,1);
    end
    
elseif numel(data_length)==1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
else
    x=x_cell;
    y=y_cell;
    
end


% ylabel を作成
ylabel(label_y);
% xlabel を作成
xlabel(label_x);
box(axes1,'on');
% 残りの座標軸プロパティの設定
set(axes1,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
% legend を作成
if num_files<20
    legend1=legend(legend_txt,'Box',legend_box,'FontName',font_name);
end
hold off
colororder(color_config)
if strcmp(fig_copy,'on')==1  % クリップボードにコピーする
    % 現在の図の背景などを透明色に指定
set(gcf,'InvertHardcopy','off')
set(gcf,'Color','none');
set(gca,'Color','none');
% クリップボードにコピーして背景などを元に戻す
try
print('-clipboard','-dmeta')
catch
    print('-clipboard','-dbitmap')
end
set(gcf,'InvertHardcopy','on')
set(gcf,'Color',[1,1,1]);
set(gca,'Color','[1,1,1]');
end
set(gcf,'visible','on')
%%  データ及び図の保存
tr=strcmp(single_save_id,'on');
tr_multi=strcmp(multiple_save_id,'on');
tr_func_back_up=strcmp(func_back_up,'on');
tr_func_back_up_multi=strcmp(multiple_func_back_up,'on');
if tr==1 || tr_multi==1 || tr_func_back_up==1 || tr_func_back_up_multi==1
    if isempty(save_dir)==1
        save_dir=oldfolder;
        cd(save_dir)
    end
    try
        cd(save_dir);
        save_dir=cd;
    catch
        cd(oldfolder)
        mkdir(save_dir);
        cd(save_dir);
        save_dir=cd;
    end
    
    
    previous_fig_file=dir(['*',fig_file_name,'*.fig']);
    previous_calc_file=dir(['*',calc_file_name,'*.mat']);
    previous_func_back_up=dir(['*',func_name,'*',]);
    
    if tr_multi==1 || tr_func_back_up_multi==1
        
        
        file_count_fig=sprintf(file_count_method,length(previous_fig_file)+1);
        file_count_calc=sprintf(file_count_method,length(previous_calc_file)+1);
        file_count_func=sprintf(file_count_method,length(previous_func_back_up)+1);
        
        fig_file_name=[fig_file_name,'_',file_count_fig];
        calc_file_name=[calc_file_name,'_',file_count_calc];
        func_back_up_file_name=[func_back_up_file_name,'_',file_count_func];
    end
    
    
    
    
    calc_file_name=append("'",calc_file_name,"'");
    
    for k=1:length(save_variables)
        
        
        save_variables_tmp1=(save_variables{k});
        
        save_variables_tmp{k}=append(",","'",save_variables{k},"'");
        expression_assign=join(["assignin('caller'",save_variables_tmp{k},",",save_variables_tmp1,")"]);
        eval(expression_assign)
        
    end
    save_variables=join(string(save_variables_tmp));
    
    expression=join(["save(",calc_file_name,save_variables,")"]);
    % expression...save( 'raw_data' , 'x','y','iv_files' )
    
    if isempty(previous_fig_file)==1
        savefig(fig_file_name);
    end
    
    if isempty(previous_calc_file)==1
        evalin('caller',expression)  % save the variales by evalin
    end
    
    
end


if tr_func_back_up==1 || tr_func_back_up_multi==1
    copyfile([func_dir,'.m'],[save_dir,'\',func_back_up_file_name,'.txt'])
end



cd(oldfolder)



    
    
end

% Raman
function [x,y,raman_files,legend_txt] = getraman_uec(varargin)
    %% control parameters
% Data files
ext='csv';     %データの拡張子
data_dir=[];   %データのディレクトリ
include_and=[];  % example include_and=[{'dark'},{'600'}] or [] ...both 600 and dark
exclude_and=[]; % example add=[{'dark'},{'600'}] or [] ...not both 600 and dark
include_or=[]; % example include_or=[{'dark'},{'600'}] or [] ..either 600 or dark
exclude_or=[];% example exclude_or=[{'dark'},{'600'}] or []..neither 600 nor dark
intervals=[]; % [] or [a,b].... Only a*x+b th data which satisfies above conditions are extracted  
move_data_dir=[];

% figures
font_name='arial';
label_y='Intensity [-]';   %y軸タイトル
label_x='Raman shift [cm^-^1]';   %x軸タイトル
font_sz=20;   %図のフォントサイズ
line_width=1.5;   %プロット線の太さ
legend_box='off';   %凡例のボックスの表示('on')、非表示('off')
fig_copy='on';    %図のクリップボードへのコピーあり('on')、なし('off')

% save
save_dir=[];    %計算結果等を保存するディレクトリ
single_save_id='off';    %セーブあり('on')、なし('off')、上書きなし
multiple_save_id='off';    %連番をつけたうえでのセーブあり('on')、なし('off')
func_back_up='off';    %使用した関数のバックアップ(txt形式)あり('on')、なし('off')
multiple_func_back_up='off';    %連番での関数のバックアップ(txt形式)あり('on')、なし('off')

fig_file_name='raw_data'; %保存する図のタイトル
calc_file_name='raw_data'; %保存するデータのタイトル
save_variables={'x','y','xrd_files'}; %保存する変数
file_count_method='%05i';  %連番の設定

% other


color_config=[        0         0    1.0000
    1.0000         0         0
    0    0.5000         0
    1.0000         0    1.0000
    0.5000         0    1.0000
    0    0.7500    1.0000
    1.0000    0.5000         0
    0.5000    0.7500         0
    1.0000    0.7500    0.7500
    0.5000    0.5000    1.0000
    0.7500    0.7500         0
    0.5000    0.5000    0.5000
    0.7500         0         0
    0         0         0];
% desing the color matrix. column corresponds to R, G, B





%% Setting
% data directly
oldfolder=cd;
if isempty(data_dir)~=1
    
    try
        cd(data_dir);
        
    catch
        mkdir(data_dir);
        cd(data_dir);
        
    end
else
    data_dir=cd;
    cd(data_dir)
end

func_dir=mfilename('fullpath');
func_name=mfilename;
previous_func_back_up=dir(['*',func_name,'*',]);
func_back_up_file_name=[func_name,'_back_up'];

num_ext=length(ext);

figure1 = figure('Color',[1 1 1],'visible','off');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
data_files_all=dir(['*',ext]);
data_files_all=data_files_all(~ismember({data_files_all.name},{previous_func_back_up.name}));
%% data search, filtering
% searching specific (and,or) data
if isempty(varargin)==1
    
    
    % and search
    % 1) include
    if isempty(include_and)==1
        idx_include_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(include_and)
            if i==1
                idx_include_and=ones(1,length(data_files_all));
            end
            include_and_tmp=['*',include_and{i},'*'];
            include_and_files=dir(include_and_tmp);
            idx_include_and_tmp=ismember({data_files_all.name},{include_and_files.name});
            idx_include_and=(idx_include_and.*idx_include_and_tmp);
            
        end
        
    end
    
    
    % 2) exclude
    if isempty(exclude_and)==1
        idx_exclude_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_and)
            if i==1
                idx_exclude_and=ones(1,length(data_files_all));
            end
            exclude_and_tmp=['*',exclude_and{i},'*'];
            exclude_and_files=dir(exclude_and_tmp);
            idx_exclude_and_tmp=ismember({data_files_all.name},{exclude_and_files.name});
            idx_exclude_and=(idx_exclude_and.*idx_exclude_and_tmp);
            
        end
        idx_exclude_and=~idx_exclude_and;
    end
    idx_and=((idx_include_and.*idx_exclude_and)>=1)';
    
    
    % or search
    % 1) include
    if isempty(include_or)==1
        idx_include_or=ones(1,length(data_files_all));
        
    else
        
        for i=1:length(include_or)
            if i==1
                idx_include_or=ones(1,length(data_files_all));
            end
            include_or_tmp=['*',include_or{i},'*'];
            include_or_files=dir(include_or_tmp);
            idx_include_or_tmp=~ismember({data_files_all.name},{include_or_files.name});
            idx_include_or=(idx_include_or.*idx_include_or_tmp);
            
        end
        idx_include_or=~idx_include_or;
    end
    
    % 2) exclude
    if isempty(exclude_or)==1
        idx_exclude_or=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_or)
            if i==1
                idx_exclude_or=logical(length(data_files_all));
            end
            exclude_or_tmp=['*',exclude_or{i},'*'];
            exclude_or_files=dir(exclude_or_tmp);
            idx_exclude_or_tmp=~ismember({data_files_all.name},{exclude_or_files.name});
            idx_exclude_or=(idx_exclude_or.*idx_exclude_or_tmp);
            
        end
    end
    
    idx_or=((idx_include_or.*idx_exclude_or)>=1)';
    
    
    idx2see=(idx_and.*idx_or)>=1;
    
    raman_files=data_files_all(idx2see);
    
    
else
    file_name=varargin;
    data_files_all(1:end)=[];
    raman_files=data_files_all;
    for j=1:length(file_name)
        
        raman_files(j)=dir(file_name{j});
    end
    
end

if isempty(intervals)~=1
    max_num_files=length(raman_files);
    
    counting_number=0;
    tmp_files=dir(ext);
    while max_num_files>=intervals(1)*counting_number+intervals(2)
        tmp_files(counting_number+1,1)=raman_files(intervals(1)*counting_number+intervals(2));
        counting_number=counting_number+1;
    end
    
    raman_files=tmp_files;
end
num_files=length(raman_files);

if num_files==0
    disp('No data satisfies the conditions you set')
    close gcf
    return
end
%%  データのインポート(以下個々のデータでコード生成して張り付け、ループにする)   
    for i=1:num_files
       
            



%IMPORTFILE2 テキスト ファイルからデータをインポート
%  UNTITLED = IMPORTFILE2(FILENAME) は既定の選択に関してテキスト ファイル FILENAME
%  からデータを読み取ります。  データを table として返します。
%
%  UNTITLED = IMPORTFILE2(FILE, DATALINES) はテキスト ファイル FILENAME
%  の指定された行区間のデータを読み取ります。DATALINES
%  を正の整数スカラーとして指定するか、行区間が不連続の場合は正の整数スカラーからなる N 行 2 列の配列として指定します。
%
%  例:
%  Untitled = importfile2("C:\Users\YuyaNagai\Desktop\Lab\Research\XRD\4mix\XRD\20200715_01_1.ras", [319, 3819]);
%
%  READTABLE も参照してください。
%
% MATLAB からの自動生成日: 2020/07/16 12:10:33

%% 入力の取り扱い

if nargin < 2
    dataLines = [20, Inf];
end

%% インポート オプションの設定およびデータのインポート
opts = delimitedTextImportOptions("NumVariables", 2);

% 範囲と区切り記号の指定
opts.DataLines = dataLines;
opts.Delimiter = ",";

% 列名と型の指定
opts.VariableNames = ["TITLE", "VarName2"];
opts.VariableTypes = ["double", "double"];

% ファイル レベルのプロパティを指定
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% データのインポート
tmp_data = readtable(raman_files(i).name, opts);



   %% 出力型への変換
% if i>1
%     data{i}(:,1)=[];
% end

nodata_ind=ismissing(tmp_data(:,2));

tmp_data(nodata_ind,:)=[];

data{i} =(table2array(tmp_data));
nan_idx = isnan(data{i});
data{i}(nan_idx,:) = [];

x_cell{i}=data{i}(:,1);
y_cell{i}=data{i}(:,2);
 data_length(i)=length(x_cell{i});

% all_data=cell2mat(data);
% x=all_data(:,1);
% y=all_data(:,2:end);

  %% グラフ描画
    plot(x_cell{i},y_cell{i},'LineWidth',1.5)
    
    % ファイル名から凡例を作成
    legend_tmp=raman_files(i).name;
    L{i}=find(legend_tmp=='_');
    legend_tmp(L{i})=' ';
    legend_tmp(end-num_ext:end)=[]; % ファイル名の拡張子部分の削除
    legend_txt{i}=legend_tmp;

 %% データ移動(コピー)
 if isempty(move_data_dir)==0
     
     try
         cd(move_data_dir)
         move_data_dir=cd;
         copyfile([raman_files(i).folder,'\',raman_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new','.',ext])
     catch
         cd(oldfolder)
         mkdir(move_data_dir);
         cd(move_data_dir);
         move_data_dir=cd;
         copyfile([raman_files(i).folder,'\',raman_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new.','.',ext])
     end
 
 end
    end
 
    %% グラフ及び出力値の調整
if length(unique(data_length))==1 && numel(data_length)>1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
    for l=1:num_files-1
        x_id(:,l)=x(:,l+1)-x(:,l);
    end
    if sum(x_id(:))==0
        x=x(:,1);
    end
    
elseif numel(data_length)==1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
else
    x=x_cell;
    y=y_cell;
    
end


% ylabel を作成
ylabel(label_y);
% xlabel を作成
xlabel(label_x);
box(axes1,'on');
% 残りの座標軸プロパティの設定
set(axes1,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
% legend を作成
if num_files<20
    legend1=legend(legend_txt,'Box',legend_box,'FontName',font_name);
end
hold off
colororder(color_config)
if strcmp(fig_copy,'on')==1  % クリップボードにコピーする
    % 現在の図の背景などを透明色に指定
set(gcf,'InvertHardcopy','off')
set(gcf,'Color','none');
set(gca,'Color','none');
% クリップボードにコピーして背景などを元に戻す
try
print('-clipboard','-dmeta')
catch
    print('-clipboard','-dbitmap')
end
set(gcf,'InvertHardcopy','on')
set(gcf,'Color',[1,1,1]);
set(gca,'Color','[1,1,1]');
end
set(gcf,'visible','on')
%%  データ及び図の保存
tr=strcmp(single_save_id,'on');
tr_multi=strcmp(multiple_save_id,'on');
tr_func_back_up=strcmp(func_back_up,'on');
tr_func_back_up_multi=strcmp(multiple_func_back_up,'on');
if tr==1 || tr_multi==1 || tr_func_back_up==1 || tr_func_back_up_multi==1
    if isempty(save_dir)==1
        save_dir=oldfolder;
        cd(save_dir)
    end
    try
        cd(save_dir);
        save_dir=cd;
    catch
        cd(oldfolder)
        mkdir(save_dir);
        cd(save_dir);
        save_dir=cd;
    end
    
    
    previous_fig_file=dir(['*',fig_file_name,'*.fig']);
    previous_calc_file=dir(['*',calc_file_name,'*.mat']);
    previous_func_back_up=dir(['*',func_name,'*',]);
    
    if tr_multi==1 || tr_func_back_up_multi==1
        
        
        file_count_fig=sprintf(file_count_method,length(previous_fig_file)+1);
        file_count_calc=sprintf(file_count_method,length(previous_calc_file)+1);
        file_count_func=sprintf(file_count_method,length(previous_func_back_up)+1);
        
        fig_file_name=[fig_file_name,'_',file_count_fig];
        calc_file_name=[calc_file_name,'_',file_count_calc];
        func_back_up_file_name=[func_back_up_file_name,'_',file_count_func];
    end
    
    
    
    
    calc_file_name=append("'",calc_file_name,"'");
    
    for k=1:length(save_variables)
        
        
        save_variables_tmp1=(save_variables{k});
        
        save_variables_tmp{k}=append(",","'",save_variables{k},"'");
        expression_assign=join(["assignin('caller'",save_variables_tmp{k},",",save_variables_tmp1,")"]);
        eval(expression_assign)
        
    end
    save_variables=join(string(save_variables_tmp));
    
    expression=join(["save(",calc_file_name,save_variables,")"]);
    % expression...save( 'raw_data' , 'x','y','iv_files' )
    
    if isempty(previous_fig_file)==1
        savefig(fig_file_name);
    end
    
    if isempty(previous_calc_file)==1
        evalin('caller',expression)  % save the variales by evalin
    end
    
    
end


if tr_func_back_up==1 || tr_func_back_up_multi==1
    copyfile([func_dir,'.m'],[save_dir,'\',func_back_up_file_name,'.txt'])
end



cd(oldfolder)



    
    
end

% PEIS
function [x,y,z,imp_files,imp_abs,freq,phase_shift] = getimp(varargin)
%% control parameters
% Data files
ext='txt';     %データの拡張子
data_dir=[];   %データのディレクトリ
include_and=[];  % example include_and=[{'dark'},{'600'}] or [] ...both 600 and dark
exclude_and=[]; % example add=[{'dark'},{'600'}] or [] ...not both 600 and dark
include_or=[]; % example include_or=[{'dark'},{'600'}] or [] ..either 600 or dark
exclude_or=[];% example exclude_or=[{'dark'},{'600'}] or []..neither 600 nor dark
intervals=[]; % [] or [a,b].... Only a*x+b th data which satisfies above conditions are extracted  
move_data_dir=[];%'C:\Users\YuyaNagai\Desktop\Lab\Research\test'

% figures
font_name='arial';
label_y="-Z'' [\Omega]";   %y軸タイトル
label_x="Z' [\Omega]";   %x軸タイトル
label_y2="Phase shift [degree]";
label_x2="Frequency [hZ]";
label_y3 = "Impedance [\Omega]";
label_x3 =label_x2;
font_sz=30;   %図のフォントサイズ
line_width=2;   %プロット線の太さ
legend_box='off';   %凡例のボックスの表示('on')、非表示('off')
fig_copy='off';    %図のクリップボードへのコピーあり('on')、なし('off')

% save
save_dir=['C:\Users\YuyaNagai\Desktop\Lab\Research\test\test'];    %計算結果等を保存するディレクトリ
single_save_id='off';    %セーブあり('on')、なし('off')、上書きなし
multiple_save_id='off';    %連番をつけたうえでのセーブあり('on')、なし('off')
func_back_up='off';    %使用した関数のバックアップ(txt形式)あり('on')、なし('off')
multiple_func_back_up='off';    %連番での関数のバックアップ(txt形式)あり('on')、なし('off')

fig_file_name='Nyquist plot'; %保存する図のタイトル
fig_file_name2='Bode plot (phase)';
fig_file_name3 = 'Bode plot (imp)';
calc_file_name='raw_data'; %保存するデータのタイトル
save_variables={'x','y','imp_files'}; %保存する変数
file_count_method='%05i';  %連番の設定

% other

% pH=13.61;

color_config=[        0         0    1.0000
    1.0000         0         0
    0    0.5000         0
    1.0000         0    1.0000
    0.5000         0    1.0000
    0    0.7500    1.0000
    1.0000    0.5000         0
    0.5000    0.7500         0
    1.0000    0.7500    0.7500
    0.5000    0.5000    1.0000
    0.7500    0.7500         0
    0.5000    0.5000    0.5000
    0.7500         0         0
    0         0         0];
% desing the color matrix. column corresponds to R, G, B





%% Setting
% data directly
oldfolder=cd;
if isempty(data_dir)~=1
    
    try
        cd(data_dir);
        
    catch
        mkdir(data_dir);
        cd(data_dir);
        
    end
else
    data_dir=cd;
    cd(data_dir)
end

func_dir=mfilename('fullpath');
func_name=mfilename;
previous_func_back_up=dir(['*',func_name,'*',]);
func_back_up_file_name=[func_name,'_back_up'];

num_ext=length(ext);

figure1 = figure('Color',[1 1 1],'visible','off');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
figure2 = figure('Color',[1 1 1],'visible','off');
axes2 = axes('Parent',figure2,'Xscale','log');
hold(axes2,'on');
figure3 = figure('Color',[1 1 1],'visible','off');
axes3 = axes('Parent',figure3,'Xscale','log','Yscale','log');
hold(axes3,'on');
data_files_all=dir(['*',ext]);
data_files_all=data_files_all(~ismember({data_files_all.name},{previous_func_back_up.name}));
%% data search, filtering
% searching specific (and,or) data
if isempty(varargin)==1
    
    
    % and search
    % 1) include
    if isempty(include_and)==1
        idx_include_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(include_and)
            if i==1
                idx_include_and=ones(1,length(data_files_all));
            end
            include_and_tmp=['*',include_and{i},'*'];
            include_and_files=dir(include_and_tmp);
            idx_include_and_tmp=ismember({data_files_all.name},{include_and_files.name});
            idx_include_and=(idx_include_and.*idx_include_and_tmp);
            
        end
        
    end
    
    
    % 2) exclude
    if isempty(exclude_and)==1
        idx_exclude_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_and)
            if i==1
                idx_exclude_and=ones(1,length(data_files_all));
            end
            exclude_and_tmp=['*',exclude_and{i},'*'];
            exclude_and_files=dir(exclude_and_tmp);
            idx_exclude_and_tmp=ismember({data_files_all.name},{exclude_and_files.name});
            idx_exclude_and=(idx_exclude_and.*idx_exclude_and_tmp);
            
        end
        idx_exclude_and=~idx_exclude_and;
    end
    idx_and=((idx_include_and.*idx_exclude_and)>=1)';
    
    
    % or search
    % 1) include
    if isempty(include_or)==1
        idx_include_or=ones(1,length(data_files_all));
        
    else
        
        for i=1:length(include_or)
            if i==1
                idx_include_or=ones(1,length(data_files_all));
            end
            include_or_tmp=['*',include_or{i},'*'];
            include_or_files=dir(include_or_tmp);
            idx_include_or_tmp=~ismember({data_files_all.name},{include_or_files.name});
            idx_include_or=(idx_include_or.*idx_include_or_tmp);
            
        end
        idx_include_or=~idx_include_or;
    end
    
    % 2) exclude
    if isempty(exclude_or)==1
        idx_exclude_or=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_or)
            if i==1
                idx_exclude_or=logical(length(data_files_all));
            end
            exclude_or_tmp=['*',exclude_or{i},'*'];
            exclude_or_files=dir(exclude_or_tmp);
            idx_exclude_or_tmp=~ismember({data_files_all.name},{exclude_or_files.name});
            idx_exclude_or=(idx_exclude_or.*idx_exclude_or_tmp);
            
        end
    end
    
    idx_or=((idx_include_or.*idx_exclude_or)>=1)';
    
    
    idx2see=(idx_and.*idx_or)>=1;
    
    imp_files=data_files_all(idx2see);
    
    
else
    file_name=varargin;
    data_files_all(1:end)=[];
    imp_files=data_files_all;
    for j=1:length(file_name)
        
        imp_files(j)=dir(file_name{j});
    end
    
end

if isempty(intervals)~=1
    max_num_files=length(imp_files);
    
    counting_number=0;
    tmp_files=dir(ext);
    while max_num_files>=intervals(1)*counting_number+intervals(2)
        tmp_files(counting_number+1,1)=imp_files(intervals(1)*counting_number+intervals(2));
        counting_number=counting_number+1;
    end
    
    imp_files=tmp_files;
end
num_files=length(imp_files);

if num_files==0
    disp('No data satisfies the conditions you set')
    close gcf
    return
end
%%  データのインポート(以下個々のデータでコード生成して張り付け、ループにする)
for i=1:num_files
 %% 変数を初期化します。
% if nargin<=2
    startRow = 18;
    endRow = inf;
% end

%% データの列をテキストとして読み取る:
% 詳細は TEXTSCAN のドキュメンテーションを参照してください。
formatSpec = '%10s%11s%12s%11s%s%[^\n\r]';

%% テキスト ファイルを開きます。
fileID = fopen(imp_files(i).name,'r');

%% データの列を書式設定に従って読み取ります。
% この呼び出しは、このコードの生成に使用されたファイルの構造に基づいています。別のファイルでエラーが発生する場合は、インポート ツールからコードの再生成を試みてください。
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% テキスト ファイルを閉じます。
fclose(fileID);

%% 数値テキストを含む列の内容を数値に変換します。
% 非数値テキストを NaN で置き換えます。
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5]
    % 入力 cell 配列のテキストを数値に変換します。非数値テキストを NaN で置き換えました。
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % 数値でない接頭辞と接尾辞を検出して削除する正規表現を作成します。
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % 桁区切り以外の場所でコンマが検出されました。
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % 数値テキストを数値に変換します。
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% 非数値セルを次の値で置き換え NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 非数値セルを検索
raw(R) = {NaN}; % 非数値セルを置き換え

%% 出力変数の作成
raw_data{i} = cell2mat(raw);
    
 
    
    
    x_cell_tmp=raw_data{i}(:,2);
    y_cell_tmp=-raw_data{i}(:,3);
    freq_tmp = raw_data{i}(:,1);
    phase_shift_tmp = raw_data{i}(:,5);
    imp_abs_tmp = raw_data{i}(:,4);
    x_cell_tmp(isnan(x_cell_tmp))=[];
    y_cell_tmp(isnan(y_cell_tmp))=[];
    freq_tmp(isnan(freq_tmp))=[];
    phase_shift_tmp(isnan(phase_shift_tmp))=[];
    imp_abs_tmp(isnan(imp_abs_tmp))=[];
    
    x_cell{i}=x_cell_tmp;
    y_cell{i}=y_cell_tmp;
    z_cell{i}=trapz(x_cell_tmp,y_cell_tmp);
    freq_cell{i} = freq_tmp;
    phase_shift_cell{i}=phase_shift_tmp;
    imp_abs_cell{i} = imp_abs_tmp;
    
    data_length(i)=length(x_cell{i});
    %% グラフ描画
    plot(axes1,x_cell{i},y_cell{i},'O-','LineWidth',line_width)
    plot(axes2,freq_cell{i},phase_shift_cell{i},'O-','LineWidth',line_width)
    plot(axes3,freq_cell{i},imp_abs_cell{i},'O-','LineWidth',line_width)
    % ファイル名から凡例を作成
    legend_tmp=imp_files(i).name;
    L{i}=find(legend_tmp=='_');
    legend_tmp(L{i})=' ';
    legend_tmp(end-num_ext:end)=[]; % ファイル名の拡張子部分の削除
    legend_txt{i}=legend_tmp;
    
 %% データ移動(コピー)
 
 if isempty(move_data_dir)==0
     
     try
         cd(move_data_dir)
         move_data_dir=cd;
         copyfile([imp_files(i).folder,'\',imp_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new','.',ext])
     catch
         cd(oldfolder)
         mkdir(move_data_dir);
         cd(move_data_dir);
         move_data_dir=cd;
         copyfile([imp_files(i).folder,'\',imp_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new.','.',ext])
     end
 end
end
%% グラフ及び出力値の調整
if length(unique(data_length))==1 && numel(data_length)>1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
    z=cell2mat(z_cell);
    freq = cell2mat(freq_cell);
    phase_shift = cell2mat(phase_shift_cell);
    imp_abs = cell2mat(imp_abs_cell);
    
    for l=1:num_files-1
        x_id(:,l)=x(:,l+1)-x(:,l);
    end
    if sum(abs(x_id(:)))==0
        x=x(:,1);
    end
    
elseif numel(data_length)==1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
    z=cell2mat(z_cell);
    freq = cell2mat(freq_cell);
    phase_shift = cell2mat(phase_shift_cell);
    imp_abs = cell2mat(imp_abs_cell);
else
    x=x_cell;
    y=y_cell;
    z=z_cell;
    freq = (freq_cell);
    phase_shift = (phase_shift_cell);
    imp_abs = (imp_abs_cell);    
end


% ylabel を作成
ylabel(axes1,label_y);
ylabel(axes2,label_y2);
ylabel(axes3,label_y3);
% xlabel を作成
xlabel(axes1,label_x);
xlabel(axes2,label_x2);
xlabel(axes3,label_x3);

box(axes1,'on');
box(axes2,'on');
box(axes3,'on');
% 残りの座標軸プロパティの設定
set(axes1,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
set(axes2,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
set(axes3,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
% legend を作成
if num_files<20
    legend1=legend(axes1,legend_txt,'Box',legend_box,'FontName',font_name);
    legend2=legend(axes2,legend_txt,'Box',legend_box,'FontName',font_name);
    legend3=legend(axes3,legend_txt,'Box',legend_box,'FontName',font_name);
end
hold(axes1, 'off');hold(axes2, 'off');hold(axes3, 'off')

colororder(figure1,color_config);colororder(figure2,color_config);colororder(figure3,color_config)

if strcmp(fig_copy,'on')==1  % クリップボードにコピーする
    % 現在の図の背景などを透明色に指定
set(figure1,'InvertHardcopy','off');set(figure2,'InvertHardcopy','off');set(figure3,'InvertHardcopy','off')
set(figure1,'Color','none');set(figure2,'Color','none');set(figure3,'Color','none');
set(axes1,'Color','none');set(axes2,'Color','none');set(axes3,'Color','none');
% クリップボードにコピーして背景などを元に戻す
try
print(figure1,'-clipboard','-dmeta');print(figure2,'-clipboard','-dmeta');print(figure3,'-clipboard','-dmeta')
catch
    print(figure1,'-clipboard','-dbitmap');print(figure2,'-clipboard','-dbitmap');print(figure3,'-clipboard','-dbitmap')
end
set(figure1,'InvertHardcopy','on');set(figure2,'InvertHardcopy','on');set(figure3,'InvertHardcopy','on');
set(figure1,'Color',[1,1,1]);set(figure2,'Color',[1,1,1]);set(figure3,'Color',[1,1,1]);
set(axes1,'Color','[1,1,1]');set(axes2,'Color','[1,1,1]');set(axes3,'Color','[1,1,1]');
end
set(figure1,'visible','on'); set(figure2,'visible','on');set(figure3,'visible','on');
%%  データ及び図の保存
tr=strcmp(single_save_id,'on');
tr_multi=strcmp(multiple_save_id,'on');
tr_func_back_up=strcmp(func_back_up,'on');
tr_func_back_up_multi=strcmp(multiple_func_back_up,'on');
if tr==1 || tr_multi==1 || tr_func_back_up==1 || tr_func_back_up_multi==1
    if isempty(save_dir)==1
        save_dir=oldfolder;
        cd(save_dir)
    end
    try
        cd(save_dir);
        save_dir=cd;
    catch
        cd(oldfolder)
        mkdir(save_dir);
        cd(save_dir);
        save_dir=cd;
    end
    
    
    previous_fig_file=dir(['*',fig_file_name,'*.fig']);
    previous_fig_file2=dir(['*',fig_file_name2,'*.fig']);
    previous_fig_file3=dir(['*',fig_file_name3,'*.fig']);
    previous_calc_file=dir(['*',calc_file_name,'*.mat']);
    previous_func_back_up=dir(['*',func_name,'*',]);
    
    if tr_multi==1 || tr_func_back_up_multi==1
        
        
        file_count_fig=sprintf(file_count_method,length(previous_fig_file)+1);
        file_count_fig2=sprintf(file_count_method,length(previous_fig_file2)+1);
        file_count_fig3=sprintf(file_count_method,length(previous_fig_file3)+1);
        file_count_calc=sprintf(file_count_method,length(previous_calc_file)+1);
        file_count_func=sprintf(file_count_method,length(previous_func_back_up)+1);
        
        fig_file_name=[fig_file_name,'_',file_count_fig];
        fig_file_name2=[fig_file_name2,'_',file_count_fig2];
        fig_file_name3=[fig_file_name3,'_',file_count_fig3];
        calc_file_name=[calc_file_name,'_',file_count_calc];
        func_back_up_file_name=[func_back_up_file_name,'_',file_count_func];
    end
    
    
    
    
    calc_file_name=append("'",calc_file_name,"'");
    
    for k=1:length(save_variables)
        
        
        save_variables_tmp1=(save_variables{k});
        
        save_variables_tmp{k}=append(",","'",save_variables{k},"'");
        expression_assign=join(["assignin('caller'",save_variables_tmp{k},",",save_variables_tmp1,")"]);
        eval(expression_assign)
        
    end
    save_variables=join(string(save_variables_tmp));
    
    expression=join(["save(",calc_file_name,save_variables,")"]);
    % expression...save( 'raw_data' , 'x','y','iv_files' )
    
    if isempty(previous_fig_file)==1
        savefig(fig_file_name);
    end
    if isempty(previous_fig_file2)==1
        savefig(fig_file_name2);
    end
    if isempty(previous_fig_file3)==1
        savefig(fig_file_name3);
    end
    if isempty(previous_calc_file)==1
        evalin('caller',expression)  % save the variales by evalin
    end
    
    
end


if tr_func_back_up==1 || tr_func_back_up_multi==1
    copyfile([func_dir,'.m'],[save_dir,'\',func_back_up_file_name,'.txt'])
end



cd(oldfolder)





end

% PEC
function [x,y,iv_files] = getiv(varargin)
%% control parameters
% Data files
ext='txt';     %データの拡張子
data_dir=[];   %データのディレクトリ
include_and=[];  % example include_and=[{'dark'},{'600'}] or [] ...both 600 and dark
exclude_and=[{'dark'}]; % example add=[{'dark'},{'600'}] or [] ...not both 600 and dark
include_or=[]; % example include_or=[{'dark'},{'600'}] or [] ..either 600 or dark
exclude_or=[];% example exclude_or=[{'dark'},{'600'}] or []..neither 600 nor dark
intervals=[]; % [] or [a,b].... Only a*x+b th data which satisfies above conditions are extracted  
move_data_dir=[];%'C:\Users\YuyaNagai\Desktop\Lab\Research\test'

% figures
font_name='arial';
label_y='Current density [mA/cm^2]';   %y軸タイトル
label_x='Potential / V vs RHE';   %x軸タイトル
font_sz=20;   %図のフォントサイズ
line_width=1.5;   %プロット線の太さ
legend_box='off';   %凡例のボックスの表示('on')、非表示('off')
fig_copy='off';    %図のクリップボードへのコピーあり('on')、なし('off')

% save
save_dir=['C:\Users\YuyaNagai\Desktop\Lab\Research\test\test'];    %計算結果等を保存するディレクトリ
single_save_id='off';    %セーブあり('on')、なし('off')、上書きなし
multiple_save_id='off';    %連番をつけたうえでのセーブあり('on')、なし('off')
func_back_up='off';    %使用した関数のバックアップ(txt形式)あり('on')、なし('off')
multiple_func_back_up='off';    %連番での関数のバックアップ(txt形式)あり('on')、なし('off')

fig_file_name='raw_data'; %保存する図のタイトル
calc_file_name='raw_data'; %保存するデータのタイトル
save_variables={'x','y','iv_files'}; %保存する変数
file_count_method='%05i';  %連番の設定

% other

pH=13.61;

color_config=[        0         0    1.0000
    1.0000         0         0
    0    0.5000         0
    1.0000         0    1.0000
    0.5000         0    1.0000
    0    0.7500    1.0000
    1.0000    0.5000         0
    0.5000    0.7500         0
    1.0000    0.7500    0.7500
    0.5000    0.5000    1.0000
    0.7500    0.7500         0
    0.5000    0.5000    0.5000
    0.7500         0         0
    0         0         0];
% desing the color matrix. column corresponds to R, G, B





%% Setting
% data directly
oldfolder=cd;
if isempty(data_dir)~=1
    
    try
        cd(data_dir);
        
    catch
        mkdir(data_dir);
        cd(data_dir);
        
    end
else
    data_dir=cd;
    cd(data_dir)
end

func_dir=mfilename('fullpath');
func_name=mfilename;
previous_func_back_up=dir(['*',func_name,'*',]);
func_back_up_file_name=[func_name,'_back_up'];

num_ext=length(ext);

figure1 = figure('Color',[1 1 1],'visible','off');
axes1 = axes('Parent',figure1);
hold(axes1,'on');
data_files_all=dir(['*',ext]);
data_files_all=data_files_all(~ismember({data_files_all.name},{previous_func_back_up.name}));
%% data search, filtering
% searching specific (and,or) data
if isempty(varargin)==1
    
    
    % and search
    % 1) include
    if isempty(include_and)==1
        idx_include_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(include_and)
            if i==1
                idx_include_and=ones(1,length(data_files_all));
            end
            include_and_tmp=['*',include_and{i},'*'];
            include_and_files=dir(include_and_tmp);
            idx_include_and_tmp=ismember({data_files_all.name},{include_and_files.name});
            idx_include_and=(idx_include_and.*idx_include_and_tmp);
            
        end
        
    end
    
    
    % 2) exclude
    if isempty(exclude_and)==1
        idx_exclude_and=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_and)
            if i==1
                idx_exclude_and=ones(1,length(data_files_all));
            end
            exclude_and_tmp=['*',exclude_and{i},'*'];
            exclude_and_files=dir(exclude_and_tmp);
            idx_exclude_and_tmp=ismember({data_files_all.name},{exclude_and_files.name});
            idx_exclude_and=(idx_exclude_and.*idx_exclude_and_tmp);
            
        end
        idx_exclude_and=~idx_exclude_and;
    end
    idx_and=((idx_include_and.*idx_exclude_and)>=1)';
    
    
    % or search
    % 1) include
    if isempty(include_or)==1
        idx_include_or=ones(1,length(data_files_all));
        
    else
        
        for i=1:length(include_or)
            if i==1
                idx_include_or=ones(1,length(data_files_all));
            end
            include_or_tmp=['*',include_or{i},'*'];
            include_or_files=dir(include_or_tmp);
            idx_include_or_tmp=~ismember({data_files_all.name},{include_or_files.name});
            idx_include_or=(idx_include_or.*idx_include_or_tmp);
            
        end
        idx_include_or=~idx_include_or;
    end
    
    % 2) exclude
    if isempty(exclude_or)==1
        idx_exclude_or=ones(1,length(data_files_all));
    else
        
        for i=1:length(exclude_or)
            if i==1
                idx_exclude_or=logical(length(data_files_all));
            end
            exclude_or_tmp=['*',exclude_or{i},'*'];
            exclude_or_files=dir(exclude_or_tmp);
            idx_exclude_or_tmp=~ismember({data_files_all.name},{exclude_or_files.name});
            idx_exclude_or=(idx_exclude_or.*idx_exclude_or_tmp);
            
        end
    end
    
    idx_or=((idx_include_or.*idx_exclude_or)>=1)';
    
    
    idx2see=(idx_and.*idx_or)>=1;
    
    iv_files=data_files_all(idx2see);
    
    
else
    file_name=varargin;
    data_files_all(1:end)=[];
    iv_files=data_files_all;
    for j=1:length(file_name)
        
        iv_files(j)=dir(file_name{j});
    end
    
end

if isempty(intervals)~=1
    max_num_files=length(iv_files);
    
    counting_number=0;
    tmp_files=dir(ext);
    while max_num_files>=intervals(1)*counting_number+intervals(2)
        tmp_files(counting_number+1,1)=iv_files(intervals(1)*counting_number+intervals(2));
        counting_number=counting_number+1;
    end
    
    iv_files=tmp_files;
end
num_files=length(iv_files);

if num_files==0
    disp('No data satisfies the conditions you set')
    close gcf
    return
end
%%  データのインポート(以下個々のデータでコード生成して張り付け、ループにする)
for i=1:num_files
    %% 入力の取り扱い
    dataLines = [20, Inf];
    %     % dataLines が指定されていない場合、既定値を定義します
    %     if nargin < 2
    %         dataLines = [20, Inf];
    %     end
    %
    %% インポート オプションの設定およびデータのインポート
    opts = delimitedTextImportOptions("NumVariables", 2);
    
    % 範囲と区切り記号の指定
    opts.DataLines = dataLines;
    opts.Delimiter = " ";
    
    % 列名と型の指定
    opts.VariableNames = ["PotentialV", "CurrentA"];
    opts.VariableTypes = ["double", "double"];
    
    % ファイル レベルのプロパティを指定
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts.ConsecutiveDelimitersRule = "join";
    opts.LeadingDelimitersRule = "ignore";
    
    % 変数プロパティを指定
    opts = setvaropts(opts, "PotentialV", "TrimNonNumeric", true);
    opts = setvaropts(opts, "PotentialV", "ThousandsSeparator", ",");
    
    % データのインポート
    mydata{i} = readtable(iv_files(i).name, opts);
    
    %% 出力型への変換
    raw_data{i} = table2array(mydata{i});
    
    x_cell_tmp=raw_data{i}(:,1)+0.059*pH+0.1976;
    y_cell_tmp=10^3*raw_data{i}(:,2);
    x_cell_tmp(isnan(x_cell_tmp))=[];
    y_cell_tmp(isnan(y_cell_tmp))=[];
    x_cell{i}=x_cell_tmp;
    y_cell{i}=y_cell_tmp;
    data_length(i)=length(x_cell{i});
    %% グラフ描画
    plot(x_cell{i},y_cell{i},'LineWidth',1.5)
    
    % ファイル名から凡例を作成
    legend_tmp=iv_files(i).name;
    L{i}=find(legend_tmp=='_');
    legend_tmp(L{i})=' ';
    legend_tmp(end-num_ext:end)=[]; % ファイル名の拡張子部分の削除
    legend_txt{i}=legend_tmp;
    
 %% データ移動(コピー)
 
 if isempty(move_data_dir)==0
     
     try
         cd(move_data_dir)
         move_data_dir=cd;
         copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new','.',ext])
     catch
         cd(oldfolder)
         mkdir(move_data_dir);
         cd(move_data_dir);
         move_data_dir=cd;
         copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir])
%          copyfile([iv_files(i).folder,'\',iv_files(i).name],[move_data_dir,'\',iv_files(i).name,'_new.','.',ext])
     end
 end
end
%% グラフ及び出力値の調整
if length(unique(data_length))==1 && numel(data_length)>1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
    for l=1:num_files-1
        x_id(:,l)=x(:,l+1)-x(:,l);
    end
    if sum(abs(x_id(:)))==0
        x=x(:,1);
    end
    
elseif numel(data_length)==1
    x=cell2mat(x_cell);
    y=cell2mat(y_cell);
else
    x=x_cell;
    y=y_cell;
    
end


% ylabel を作成
ylabel(label_y);
% xlabel を作成
xlabel(label_x);
box(axes1,'on');
% 残りの座標軸プロパティの設定
set(axes1,'FontSize',font_sz,'LineWidth',line_width,'FontName',font_name);
% legend を作成
if num_files<20
    legend1=legend(legend_txt,'Box',legend_box,'FontName',font_name);
end
hold off
colororder(color_config)
if strcmp(fig_copy,'on')==1  % クリップボードにコピーする
    % 現在の図の背景などを透明色に指定
set(gcf,'InvertHardcopy','off')
set(gcf,'Color','none');
set(gca,'Color','none');
% クリップボードにコピーして背景などを元に戻す
try
print('-clipboard','-dmeta')
catch
    print('-clipboard','-dbitmap')
end
set(gcf,'InvertHardcopy','on')
set(gcf,'Color',[1,1,1]);
set(gca,'Color','[1,1,1]');
end
set(gcf,'visible','on')
%%  データ及び図の保存
tr=strcmp(single_save_id,'on');
tr_multi=strcmp(multiple_save_id,'on');
tr_func_back_up=strcmp(func_back_up,'on');
tr_func_back_up_multi=strcmp(multiple_func_back_up,'on');
if tr==1 || tr_multi==1 || tr_func_back_up==1 || tr_func_back_up_multi==1
    if isempty(save_dir)==1
        save_dir=oldfolder;
        cd(save_dir)
    end
    try
        cd(save_dir);
        save_dir=cd;
    catch
        cd(oldfolder)
        mkdir(save_dir);
        cd(save_dir);
        save_dir=cd;
    end
    
    
    previous_fig_file=dir(['*',fig_file_name,'*.fig']);
    previous_calc_file=dir(['*',calc_file_name,'*.mat']);
    previous_func_back_up=dir(['*',func_name,'*',]);
    
    if tr_multi==1 || tr_func_back_up_multi==1
        
        
        file_count_fig=sprintf(file_count_method,length(previous_fig_file)+1);
        file_count_calc=sprintf(file_count_method,length(previous_calc_file)+1);
        file_count_func=sprintf(file_count_method,length(previous_func_back_up)+1);
        
        fig_file_name=[fig_file_name,'_',file_count_fig];
        calc_file_name=[calc_file_name,'_',file_count_calc];
        func_back_up_file_name=[func_back_up_file_name,'_',file_count_func];
    end
    
    
    
    
    calc_file_name=append("'",calc_file_name,"'");
    
    for k=1:length(save_variables)
        
        
        save_variables_tmp1=(save_variables{k});
        
        save_variables_tmp{k}=append(",","'",save_variables{k},"'");
        expression_assign=join(["assignin('caller'",save_variables_tmp{k},",",save_variables_tmp1,")"]);
        eval(expression_assign)
        
    end
    save_variables=join(string(save_variables_tmp));
    
    expression=join(["save(",calc_file_name,save_variables,")"]);
    % expression...save( 'raw_data' , 'x','y','iv_files' )
    
    if isempty(previous_fig_file)==1
        savefig(fig_file_name);
    end
    
    if isempty(previous_calc_file)==1
        evalin('caller',expression)  % save the variales by evalin
    end
    
    
end


if tr_func_back_up==1 || tr_func_back_up_multi==1
    copyfile([func_dir,'.m'],[save_dir,'\',func_back_up_file_name,'.txt'])
end



cd(oldfolder)





end


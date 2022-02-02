function [features_tbl,features_tbl_all,calc_para,All_data_bare_hem]=feature_ext_tn(All_data_bare_hem)

%% UV-vis

%---------- convolution---------- 
% there are some unexpected distortion or split of each peak, To solve this
% convolution whose kernel is gaussian is applied.
% Parameter below need to be decided intractively.
% Example of convolution parameter is shown below.
conv_para.func_size=100; % size of kernel function
conv_para.sigma=2.5; % width of kernel
conv_para.center=50; % center of kernel 
conv_para.graph='off'; % if on, graph of kernel is presented.
% func_kernel = 1/sqrt(2*pi*conv_para.sigma^2)*(exp(-(t - conv_para.center).^2/(2*conv_para.sigma^2)));

[processed_data_abs,calc_para.uv_vis.conv] = conv_datam(All_data_bare_hem.UV_vis.absorbance,conv_para);

absfig(All_data_bare_hem.UV_vis.wavelength,processed_data_abs)
% Peak feature selection
% 2 dominant peak intensity and positions are extracted 
R1 = 300<=All_data_bare_hem.UV_vis.wavelength & All_data_bare_hem.UV_vis.wavelength<=500;
peak_reg_x_uv_vis = All_data_bare_hem.UV_vis.wavelength(R1);
peak_reg_y_uv_vis = processed_data_abs(R1,:);

 for i=1:length(All_data_bare_hem.UV_vis.files)
[pks_uv_vis{i},locs_uv_vis{i}]=findpeaks(peak_reg_y_uv_vis(:,i),'MinPeakProminence',.005,'MinPeakDistance',100);

 end
% 
 pks_uv_vis=cell2mat(pks_uv_vis);
locs_uv_vis=cell2mat(locs_uv_vis);

 features.UV_vis_R1_pks = (pks_uv_vis(1,:))' ; features.UV_vis_R1_locs = peak_reg_x_uv_vis(locs_uv_vis(1,:));
 features.UV_vis_R2_pks = (pks_uv_vis(2,:))' ; features.UV_vis_R2_locs = peak_reg_x_uv_vis(locs_uv_vis(2,:));

 % Baseline feature selection
 R2 = 700<=All_data_bare_hem.UV_vis.wavelength & All_data_bare_hem.UV_vis.wavelength<=800;
 baseline_reg_x_UV_vis = (All_data_bare_hem.UV_vis.wavelength(R2));
 base_reg_y_UV_vis = processed_data_abs(R2,:);
 
 features.UV_vis_R3_ave_tail = (mean(base_reg_y_UV_vis))';
All_data_bare_hem.UV_vis.processed_abs = processed_data_abs;
 
%% Raman
 %Data analysis for UEC raman
 clearvars -except features All_data_bare_hem calc_para 
 
x_all = All_data_bare_hem.Raman.raman_shift;
y_all = All_data_bare_hem.Raman.intensity;
% close all
% [x,y,raman_files]=getraman;
% load('C:\Users\YuyaNagai\Desktop\Lab\Research\experimental design\Fe2O3\Ramn_UEC\1st\600-12\raw_data.mat')
% L=x>1200;
% x(L)=[];
% y(L,:)=[];
% ramanfig(x,y)
% yft=fft(y);
% 
% yft(150:end,:)=0; % low pass filter 
% y=ifft(yft,'symmetric');
% ramanfig(x,y)

% clear
% load('C:\Users\YuyaNagai\Desktop\Lab\Research\experimental design\Fe2O3\Raman_UEC\1st\600-12\raw_data.mat')

% L=x<150 | x>800;
clear R1

for i = 1:size(y_all,2)
    x = x_all(:,i);
    y = y_all(:,i);
    R1=x>1200;
x(R1)=[];
y(R1)=[];
raw_data=y;
[back_corr_y]=msbackadj(x,y,'WindowSize',25,'StepSize',40,'RegressionMethod','linear','EstimationMethod', 'quantile', 'QuantileValue', 0.001, 'SmoothMethod', 'rloess');
% [back_corr_y]=msbackadj(x,y,'WindowSize',20,'StepSize',10,'ShowPlot', 20,'RegressionMethod','linear','EstimationMethod', 'quantile', 'QuantileValue', 0.001, 'SmoothMethod', 'rloess');
L_b=logical(sum(isnan(back_corr_y),2));

back_corr_y(L_b,:)=[];
x(L_b)=[];
    


% ramanfig(x,back_corr_y)
%----- sharp peak removal by spline fit(inverse backgrond correction)------
[back_corr_y2]=msbackadj(x,back_corr_y,'WindowSize',2,'StepSize',2,'RegressionMethod','linear','EstimationMethod', 'quantile', 'QuantileValue', 0.001, 'SmoothMethod', 'rloess');
L_b2=logical(sum(isnan(back_corr_y2),2));

back_corr_y2(L_b2,:)=[];
x(L_b2)=[];
back_corr_y(L_b2,:)=[];
% ramanfig(x,back_corr_y-back_corr_y2)
L3 = x<460;

%-------- fft filtering (low pass filter)-------
% L_filt =x>450;
% y_filt = back_corr_y(L_filt,:);
% yft = (fft(y_filt));
 yft = (fft(back_corr_y));
% figure;plot(real(yft));
thres =150;
yft(thres:end,:)=0;
% yft(thres:end-thres,:)=0;
% yft(end-thres+1:end,:)=0;
defilt_y = ifft(yft,'symmetric');
% y2 = back_corr_y;
% y2(L_filt,:)=defilt_y;

% defilt_y = ifft(yft,'symmetric');
y2=defilt_y;
back_corr_y_defilt = back_corr_y-back_corr_y2;
back_corr_y_defilt(L3,:) = y2(L3,:);
% ramanfig(x,back_corr_y_defilt)



%-------- convolution------- 
% there are some unexpected distortion or split of each peak, To solve this
% convolution whose kernel is gaussian is applied.
% Parameter below need to be decided intractively.
% Example of convolution parameter is shown below.
conv_para.func_size=100; % size of kernel function
conv_para.sigma=1; % width of kernel
conv_para.center=50; % center of kernel 
conv_para.graph='off'; % if on, graph of kernel is presented.
% func_kernel = 1/sqrt(2*pi*conv_para.sigma^2)*(exp(-(t - conv_para.center).^2/(2*conv_para.sigma^2)));

[conv_data,calc_para.raman.conv] = conv_datam(back_corr_y_defilt,conv_para);
processed_data{i}=conv_data;
% processed_data(L_conv_edge_adj,:)=[];
% ramanfig(x,processed_data)
raman_shift{i} =x;
end

All_data_bare_hem.Raman.processed.intensity = processed_data;
All_data_bare_hem.Raman.processed.raman_shift = raman_shift;
% %% fft fliter low pass
% % yft = (fft(processed_data));
% % figure;plot(real(yft));
% % 
% % yft(1:30,:)=0;
% % defilt_y = ifft(yft,'symmetric');
% % % y2 = back_corr_y;
% % % y2(L_filt,:)=defilt_y;
% % 
% % % defilt_y = ifft(yft,'symmetric');
% % y2=defilt_y;
% % ramanfig(x,y2)


%----- peak section------
% Auto peak region separation
% mean_spec_processing = mean(processed_data,2);
%  Int_threshold4peak_region = 5;
%  Dist_threshold4peak_region = 10;
%  peak_index = mean_spec_processing>Int_threshold4peak_region;
%  peak_index_re = [0;peak_index;0];
%  peak_edge_idx = (diff(peak_index_re,2))==-1;
%  
%  peak_edge = x(peak_edge_idx);
% D = (reshape(peak_edge,[2,sum(peak_edge_idx)/2]))';
% peak_dist = abs(D(:,2)-D(:,1));
% D((peak_dist<Dist_threshold4peak_region),:)=[];

D=[
   
220   235     0
235   255     0
280   310     0
395   425     0
480   525     0
595   635     0
640   680     0
790   850     0
1020   1080     0
1080   1140     0 ];

% [xc,yc,zc,x_d,y_d,z_d,k]=dividespec2(x,processed_data,D);

for i =1:length(processed_data)

    for j = 1:size(D,1)
        D_min = D(j,1);
        D_max = D(j,2);
        idx_reg = (D_min <= raman_shift{i}) & (raman_shift{i}<= D_max);
        x_reg = raman_shift{i}(idx_reg);
        y_reg = processed_data{i}(idx_reg);
        z_reg = abs(trapz(x_reg,y_reg));
        [pks_tmp,idx_tmp] = max(y_reg);
        locs_tmp = x_reg(idx_tmp);

        pks_raman{j}(i,1) = pks_tmp;
        locs_raman{j}(i,1) = locs_tmp;
        z_d(i,j) = z_reg;
    end

end


% z_d=z_d';
% for j=1:length(xc)
%         [pks_raman_tmp,locs_ori{j,1}]=max(yc{j},[],1);
% %         for k=1:length(locs_ori{j})
% %           locs  
% %         locs{j}=xc{j}
% pks_raman{j,1} = pks_raman_tmp';
%    locs_raman{j,1} = (xc{j}(locs_ori{j,1}));
%    
% end    

features.raman_R1_pks=pks_raman{1};features.raman_R1_locs=locs_raman{1};features.raman_R1_area=z_d(:,1);features.raman_R1_width=z_d(:,1)./pks_raman{1};
features.raman_R2_pks=pks_raman{2};features.raman_R2_locs=locs_raman{2};features.raman_R2_area=z_d(:,2);features.raman_R2_width=z_d(:,2)./pks_raman{2};
features.raman_R3_pks=pks_raman{3};features.raman_R3_locs=locs_raman{3};features.raman_R3_area=z_d(:,3);features.raman_R3_width=z_d(:,3)./pks_raman{3};
features.raman_R4_pks=pks_raman{4};features.raman_R4_locs=locs_raman{4};features.raman_R4_area=z_d(:,4);features.raman_R4_width=z_d(:,4)./pks_raman{4};
features.raman_R5_pks=pks_raman{5};features.raman_R5_locs=locs_raman{5};features.raman_R5_area=z_d(:,5);features.raman_R5_width=z_d(:,5)./pks_raman{5};
features.raman_R6_pks=pks_raman{6};features.raman_R6_locs=locs_raman{6};features.raman_R6_area=z_d(:,6);features.raman_R6_width=z_d(:,6)./pks_raman{6};
features.raman_R7_pks=pks_raman{7};features.raman_R7_locs=locs_raman{7};features.raman_R7_area=z_d(:,7);features.raman_R7_width=z_d(:,7)./pks_raman{7};
features.raman_R8_pks=pks_raman{8};features.raman_R8_locs=locs_raman{8};features.raman_R8_area=z_d(:,8);features.raman_R8_width=z_d(:,8)./pks_raman{8};
features.raman_R9_pks=pks_raman{9};features.raman_R9_locs=locs_raman{9};features.raman_R9_area=z_d(:,9);features.raman_R9_width=z_d(:,9)./pks_raman{9};
features.raman_R10_pks=pks_raman{10};features.raman_R10_locs=locs_raman{10};features.raman_R10_area=z_d(:,10);features.raman_R10_width=z_d(:,10)./pks_raman{10};



% features.pks=cell2mat(pks_raman);
% features.locs=cell2mat(locs_raman);
% features.area = z_d;
% features.width = features.area./features.pks;
% cos_s=cossimilarity2((pks-mean(pks,2))');
 %% PEIS
 clearvars -except features All_data_bare_hem calc_para 
 
 freq_all = All_data_bare_hem.PEIS.freq;
 z_re_all = All_data_bare_hem.PEIS.z_re;
 z_im_all = All_data_bare_hem.PEIS.z_im;
imp_abs_all = All_data_bare_hem.PEIS.imp;
 phase_shift_all = All_data_bare_hem.PEIS.phase_shift;
 %--------- Data size adjustment ---------

for i=1:length(freq_all)
    min_freq_values(i) = min(freq_all{i});
   
end

min_freq = max(min_freq_values);

for i=1:length(freq_all)
    tmp_freq = freq_all{i};
    tmp_z_re = z_re_all{i};
    tmp_z_im = z_im_all{i};
    tmp_imp_abs = imp_abs_all{i};
    tmp_phase = phase_shift_all{i};
    
    freq_cell{i} = tmp_freq(tmp_freq>=min_freq);
    z_re_cell{i} = tmp_z_re(tmp_freq>=min_freq);
    z_im_cell{i} = tmp_z_im(tmp_freq>=min_freq);
    imp_abs_cell{i} = tmp_imp_abs(tmp_freq>=min_freq);
    phase_shift_cell{i} = tmp_phase(tmp_freq>=min_freq);
    data_length(i) = length(freq_cell{i});
end

if length(unique(data_length))==1 && numel(data_length)>1
    z_re=cell2mat(z_re_cell);
    z_im=cell2mat(z_im_cell);
    
    freq = cell2mat(freq_cell);
    phase_shift = cell2mat(phase_shift_cell);
    imp_abs = cell2mat(imp_abs_cell);
    
    for l=1:length(freq_all)-1
        freq_id(:,l)=freq(:,l+1)-freq(:,l);
    end
    if sum(freq_id(:))==0
        freq=freq(:,1);
    end
    
elseif numel(data_length)==1
    z_re=cell2mat(z_re_cell);
    z_im=cell2mat(z_im_cell);
    
    freq = cell2mat(freq_cell);
    phase_shift = cell2mat(phase_shift_cell);
    imp_abs = cell2mat(imp_abs_cell);
end

freq_log = log10(freq);

%------- Feature extractions----------

% phase shift data
% separate the data into three frequency region;10^-1~10^0; 10^0~10^2;10^2~

D = [   -1.09308341467501         0
              0    2.0000
         2.0000    3.9853
         ];
     

     
     for i=1:size(D,1)
                  
         Lmin=D(i,1);
         Lmax=D(i,2);
         freq_reg{i} = freq_log(Lmin<=freq_log(:,i) & freq_log(:,i)<=Lmax,:);
         
         phase_shift_reg{i} = phase_shift(Lmin<=freq_log(:,i) & freq_log(:,i)<=Lmax,:);
         [min_phase_reg_tmp,min_index_phase_reg{i}] = min(phase_shift_reg{i}); % minimum of phase shift
         min_phase_reg{i} = (min_phase_reg_tmp)';
         min_freq_shift_reg_tmp = freq_reg{i}(min_index_phase_reg{i}); % index of minimum of phase shift
         min_freq_shift_reg{i} = (min_freq_shift_reg_tmp)';
         [max_phase_reg_tmp,max_index_phase_reg{i}] = max(phase_shift_reg{i});
         max_phase_reg{i} = (max_phase_reg_tmp)';
         max_freq_shift_reg_tmp = freq_reg{i}(max_index_phase_reg{i});
          max_freq_shift_reg{i} =  (max_freq_shift_reg_tmp)';
         
         imp_abs_reg{i} = imp_abs(Lmin<=freq_log(:,i) & freq_log(:,i)<=Lmax,:);
         [min_imp_reg_tmp,min_index_imp_reg{i}] = min(imp_abs_reg{i}); % minimum of phase shift
         min_imp_reg{i} = (min_imp_reg_tmp)';
         min_freq_imp_reg_tmp = freq_reg{i}(min_index_imp_reg{i}); % index of minimum of phase shift
          min_freq_imp_reg{i} =  (min_freq_imp_reg_tmp)';
         [max_imp_reg_tmp,max_index_imp_reg{i}] = max(imp_abs_reg{i});
         max_imp_reg{i} = (max_imp_reg_tmp)';       
         max_freq_imp_reg_tmp = freq_reg{i}(max_index_imp_reg{i});
         max_freq_imp_reg{i} = (max_freq_imp_reg_tmp)';
%                   
         
         
         for j =1:length(freq_all)
             Area_phase_reg{i}(:,j) = trapz(freq_reg{i}(:,j),phase_shift_reg{i}(:,j));
             Area_imp_reg{i}(:,j) = trapz(freq_reg{i}(:,j),imp_abs_reg{i}(:,j));
         end
         
         
         
     end
     % -----  phase shift-------
features.PEIS_phase_R1_min =  min_phase_reg{1}; features.PEIS_phase_R1_min_freq =  min_freq_shift_reg{1}; % min and loc
features.PEIS_phase_R1_max =  max_phase_reg{1}; features.PEIS_phase_R1_max_freq =  max_freq_shift_reg{1}; % max and loc
features.PEIS_phase_R1_area =  (Area_phase_reg{1})'; % area
features.PEIS_phase_R2_min =  min_phase_reg{2}; features.PEIS_phase_R2_min_freq =  min_freq_shift_reg{2}; % min and loc
features.PEIS_phase_R2_max =  max_phase_reg{2}; features.PEIS_phase_R2_max_freq =  max_freq_shift_reg{2}; % max and loc
features.PEIS_phase_R2_area =  (Area_phase_reg{2})'; % area
features.PEIS_phase_R3_min =  min_phase_reg{3}; features.PEIS_phase_R3_min_freq =  min_freq_shift_reg{3}; % min and loc
features.PEIS_phase_R3_max =  max_phase_reg{3}; features.PEIS_phase_R3_max_freq =  max_freq_shift_reg{3}; % max and loc
features.PEIS_phase_R3_area =  (Area_phase_reg{3})';  %area

%-------- impedance---------
features.PEIS_imp_R1_min =  min_imp_reg{1}; features.PEIS_imp_R1_min_freq =  min_freq_imp_reg{1}; % min and loc
features.PEIS_imp_R1_max =  max_imp_reg{1}; features.PEIS_imp_R1_max_freq =  max_freq_imp_reg{1}; % max and loc
features.PEIS_imp_R1_area =  (Area_imp_reg{1})'; % area
features.PEIS_imp_R2_min =  min_imp_reg{2}; features.PEIS_imp_R2_min_freq =  min_freq_imp_reg{2}; % min and loc
features.PEIS_imp_R2_max =  max_imp_reg{2}; features.PEIS_imp_R2_max_freq =  max_freq_imp_reg{2}; % max and loc
features.PEIS_imp_R2_area =  (Area_imp_reg{2})'; % area
features.PEIS_imp_R3_min =  min_imp_reg{3}; features.PEIS_imp_R3_min_freq =  min_freq_imp_reg{3}; % min and loc
features.PEIS_imp_R3_max =  max_imp_reg{3}; features.PEIS_imp_R3_max_freq =  max_freq_imp_reg{3}; % max and loc
features.PEIS_imp_R3_area =  (Area_imp_reg{3})'; % area

 


% min_phase_reg = cell2mat(min_phase_reg');
% min_freq_shift_reg =cell2mat(min_freq_shift_reg');
% max_phase_reg = cell2mat(max_phase_reg');
% max_freq_shift_reg = cell2mat(max_freq_shift_reg');
% Area_phase_reg=cell2mat(Area_phase_reg');
% 
% min_imp_reg = cell2mat(min_imp_reg');
% min_freq_imp_reg =cell2mat(min_freq_imp_reg');
% max_imp_reg = cell2mat(max_imp_reg');
% max_freq_imp_reg = cell2mat(max_freq_imp_reg');
% Area_imp_reg=cell2mat(Area_imp_reg');
% 
% Features_phase = abs(vertcat(min_phase_reg,min_freq_shift_reg,max_phase_reg,max_freq_shift_reg,Area_phase_reg));
% Features_imp = abs(vertcat(min_imp_reg,min_freq_imp_reg,max_imp_reg,max_freq_imp_reg,Area_imp_reg));
% 
% 
% features.imp_data = (vertcat(Features_phase,Features_imp))';
%% XRD
 clearvars -except features All_data_bare_hem calc_para 
 x = All_data_bare_hem.XRD.diff_angle;
 y = All_data_bare_hem.XRD.intensity;
L=x<25;
x(L)=[];
y(L,:)=[];
 [back_corr_y]=msbackadj(x,y,'WindowSize',3,'StepSize',3,'ShowPlot', 1,.....
 'RegressionMethod','spline','EstimationMethod', 'quantile', 'QuantileValue', 0.35,'SmoothMethod', 'loess');
xrdfig(x,back_corr_y)
% [aligned_data,~,shift]=icoshift('average',back_corr_y',10);

%--- convolution --------
% there are some unexpected distortion or split of each peak, To solve this
% convolution whose kernel is gaussian is applied.
% Parameter below need to be decided intractively.
% Example of convolution parameter is shown below.
conv_para.func_size=100; % size of kernel function
conv_para.sigma=5; % width of kernel
conv_para.center=50; % center of kernel 
conv_para.graph='off'; % if on, graph of kernel is presented.
% func_kernel = 1/sqrt(2*pi*conv_para.sigma^2)*(exp(-(t - conv_para.center).^2/(2*conv_para.sigma^2)));

[conv_data,calc_para.XRD.conv] = conv_datam(back_corr_y,conv_para);
% xrdfig(x,conv_data)

%------- Peak alignment-------

x_range=[51.5,53];
[~,~,~,x_peak,y_peak]=dividespec2(x,conv_data,x_range);
[~,peak_Ind]=max(y_peak);

shift=round(mean(peak_Ind))-peak_Ind;

for i=1:length(shift)
    y_tmp=conv_data(:,i);
    if shift(i)>0
        y_tmp(end-shift(i)+1:end) = [];
        conv_data(:,i) = vertcat(zeros(shift(i),1),y_tmp);
    elseif shift(i)<0
        y_tmp(1:abs(shift(i))) = [];
        conv_data(:,i) = vertcat(y_tmp,zeros(abs(shift(i)),1));
    end
end
xrdfig(x,conv_data)
        
processed_data = conv_data;    
All_data_bare_hem.XRD.processed.intensity = processed_data;
All_data_bare_hem.XRD.processed.d_angle = x;


%-------- peak section---------
D=[
   26.5800   28.1200         0
   33.8000   35.2400         0
   35.7600   36.9400         0
   37.8800   40.0000         0
   42.7400   45.0600         0
   51.6000   53.0000         0
   54.7400   55.8600         0
   61.7600   62.9400         0
   64.2200   65.5000         0
   65.5000   67.3400         0
   70.9600   72.3400         0
   77.6600   79.9000         0];

[xc,yc,zc,x_d,y_d,z_d,k]=dividespec2(x,processed_data,D);
% xrdfig(x_d,y_d)
   
for j=1:length(xc)
        [pks_tmp,locs_ori{j,1}]=max(yc{j},[],1);
        pks_xrd{j,1} = (pks_tmp)';
%         for k=1:length(locs_ori{j})
%           locs  
%         locs{j}=xc{j}
   locs_xrd{j,1} = (xc{j}(locs_ori{j,1}));
   
end

features.xrd_R1_pks=pks_xrd{1};features.xrd_R1_locs=locs_xrd{1};
features.xrd_R2_pks=pks_xrd{2};features.xrd_R2_locs=locs_xrd{2};
features.xrd_R3_pks=pks_xrd{3};features.xrd_R3_locs=locs_xrd{3};
features.xrd_R4_pks=pks_xrd{4};features.xrd_R4_locs=locs_xrd{4};
features.xrd_R5_pks=pks_xrd{5};features.xrd_R5_locs=locs_xrd{5};
features.xrd_R6_pks=pks_xrd{6};features.xrd_R6_locs=locs_xrd{6};
features.xrd_R7_pks=pks_xrd{7};features.xrd_R7_locs=locs_xrd{7};
features.xrd_R8_pks=pks_xrd{8};features.xrd_R8_locs=locs_xrd{8};
features.xrd_R9_pks=pks_xrd{9};features.xrd_R9_locs=locs_xrd{9};
features.xrd_R10_pks=pks_xrd{10};features.xrd_R10_locs=locs_xrd{10};
features.xrd_R11_pks=pks_xrd{11};features.xrd_R11_locs=locs_xrd{11};
features.xrd_R12_pks=pks_xrd{12};features.xrd_R12_locs=locs_xrd{12};



 %% PEC
 vol = All_data_bare_hem.PEC.vol;
 current_norm = All_data_bare_hem.PEC.current_density;
 [~,index1]=min(abs(vol-1.23));

 for i =1:size(index1,2) 
pec_target(i,1) = current_norm(index1(i),i);  % cuurent @ 1.23 V
 end
 features_tbl_ori = [struct2table(features),array2table(pec_target)];
 
 features_tbl_ori.Properties.VariableNames{end} = 'PEC';
 
 %% delete small variance 
 features_tbl_all = features_tbl_ori;
 rsd =abs(std(features_tbl_all{:,:})./mean(features_tbl_all{:,:}));
 ind_rsd = abs(rsd)<10^-2;
 features_tbl_ori(:,ind_rsd) = [];
%  
%  
%  %% auto scaling 
%  scl_features_tbl = array2table((features_tbl{:,:}-mean(features_tbl{:,:}))./std(features_tbl{:,:}),'VariableNames',features_tbl.Properties.VariableNames);
 features_tbl_inv = array2table(1./features_tbl_ori{:,1:end-1});
 for i =1:size(features_tbl_inv,2)
    features_tbl_inv.Properties.VariableNames{i} = ['inv_', features_tbl_ori.Properties.VariableNames{i}];
 end

 features_tbl = [features_tbl_inv,features_tbl_ori];


end

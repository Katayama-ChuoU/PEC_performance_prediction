%% myTbl_rev1, uv/vis position changed to inverse
%% myTbl_rev2, PEIS data removed
%% myTbl_rev3, only PEIS data
%% myTbl_rev4, high correlation descriptors without PEIS

%% Data & parameter import
load('features.mat') 
load('ipt_prm.mat')
load('ML_para.mat')


%% autoscaling of the table
myTbl = features_tbl;
myTbl_matrix = myTbl.Variables;
myTbl_size = size(myTbl_matrix);
myTbl_name = myTbl.Properties.VariableNames;

for i = 1:myTbl_size(2)
    scl_myTbl_matrix(:,i) = (myTbl_matrix(:,i) - mean(myTbl_matrix(:,i)))/std(myTbl_matrix(:,i));
end

scl_myTbl = table;
scl_myTbl.Variables = scl_myTbl_matrix;
scl_myTbl.Properties.VariableNames = myTbl_name;

%% Remove the zero/low relative standard deviation
rsd = std(features_tbl{:,:})./mean(features_tbl{:,:});
rstd_thrs = ML_para.rstd_threshold;
constant_dcp_idx = abs(rsd) < abs(rstd_thrs);
scl_myTbl(:,constant_dcp_idx) = [];
scl_myTbl_matrix(:,constant_dcp_idx) = [];

%% pca
% pca and recovery
% [coeff,score,latent] = pca(mytable);
% % this is recovered, but centered
% tmp = score*coeff';
% % by adding the mean of each variable, the data is recovered.
% tmp2 = tmp + mean(mytable(:,1:4));
% % latent includes the dispersion of 1st, 2nd, 3rd, ... component. 
% % taking square root provides the rough contribution for each pc

[coeff, score, latent] = pca(scl_myTbl_matrix);
% recovered cause the data is already autoscaled.
% tmp = score * coeff';
for i = 1:size(latent)
    latent_accum(i) = sum(sqrt(latent(1:i)));
end
latent_accum = latent_accum / sum(sqrt(latent(:)));
figure(1)
plot(latent_accum, 'o')
xlabel('Number of PCA')
% title('Cumulative frequency of PCA')

num_pca = 10;
for i = 1:myTbl_size(2)
   coeff_ave(i) = sum(coeff(i,1:num_pca)' .* sqrt(latent(1:num_pca))); 
end
%coeff_ave = mean(coeff(:,1:num_pca), 2);
figure(2)
plot(coeff_ave, 'o')
xlabel('Variables number')
% title('Impacet of each variables')


% import mlreportgen.ppt.*
% ppt = Presentation('AutomatedPresentation.pptx','myTemplate2.pptx');
% masters = getMasterNames(ppt);
% layoutnames = getLayoutNames(ppt,masters{1});
% newslide = add(ppt,layoutnames{8});
% 
% figs.fig1 = figure(1);
% % figs.fig1.Units = 'centimeter';
% % figs.fig1.PaperSize = [11.99 11.29];  % Size is same as 1 slide in ppt.
% % figs.fig1.Position(3) = 11.99; % Figure size(width) is defined.
% % figs.fig1.Position(4) = 11.29;
% adjfig
% 
% % ax = gca;
% % % ax.PositionConstraint = 'outerposition';
% % % ax.Units = 'centimeter';
% % % ax.OuterPosition(3) = figs.fig1.OuterPosition(3);
% % % ax.OuterPosition(4) = figs.fig1.OuterPosition(4);
% % outerpos = ax.OuterPosition;
% % 
% % ti = ax.TightInset; 
% % left = outerpos(1) + ti(1);
% % bottom = outerpos(2) + ti(2);
% % ax_width = outerpos(3) - ti(1) - ti(3);
% % ax_height = outerpos(4) - ti(2) - ti(4);
% % ax.Position = [left bottom ax_width ax_height];
% % ax.Units = 'centimeter';
% % figs.fig1.PaperSize = [figs.fig1.PaperPosition(3) figs.fig1.PaperPosition(4)];
% saveas(figs.fig1,'fig1.fig')
% exportgraphics(figs.fig1,'fig1.emf','BackgroundColor','none')
% picture1 = Picture(('fig1.emf'));
% 
% % picture1.Width = [num2str(figs.fig1.Position(3)),'cm'];
% % picture1.Height = [num2str(figs.fig1.Position(4)),'cm'];
% picture1.Name = 'Cumulative frequency of PCA';
% 
% % add(newslide,picture1);     % Add the figure in ppt
% replace(newslide.Children(14),picture1.Name)
% % replace(newslide.Children(9),Picture("fig1.emf"))
% replace(newslide.Children(13),picture1)
% 
% % link for 1st fig
% fig_file1 = dir('fig1.fig');
% fileloc1 = [fig_file1.folder,'\',fig_file1.name];
% % fileloc = 'figs.mat';
% link1 = ExternalLink(fileloc1,'figlink1');
% p1 = Paragraph('');
% append(p1,link1);
% replace(ppt.Children(1),newslide.Children(5).Name,p1);
% 
% figs.fig2 = figure(2);
% adjfig
% saveas(figs.fig2,'fig2.fig')
% exportgraphics(figs.fig2,'fig2.emf','BackgroundColor','none')
% 
% picture2 = Picture(('fig2.emf'));
% picture2.Name = 'Impacet of each variables';
% replace(newslide.Children(16),picture1.Name)
% replace(newslide.Children(15),picture2)
% close(ppt)

%% PLS

k_fold = ipt_prm.k_fold;
nm_pls = 5;
nm_comp = 10;

X = scl_myTbl_matrix(:, 1:(end-1));
y = scl_myTbl_matrix(:, end);

[Xloadings, Yloadings, Xscores, Yscores, betaPLS, PLSPctVar] = plsregress(X, y, 10);
yfitPLS = [ones(size(X,1), 1) X] * betaPLS;

% statistical evaluation
TSS = sum((y - mean(y)).^2);
RSS_PLS = sum((y - yfitPLS).^2);
rsquaredPLS = 1 - RSS_PLS/TSS

% cross validation
% pctVar 1st row: percentage of variance expained in x, 2nd row: percentage
% of variance explained in y
% PLSmsep: 1st row : mse for x, 2nd row: mse for y with component n
[Xl, Yl, Xs, Ys, beta, pctVar, PLSmsep] = plsregress(X, y, 10, 'CV', k_fold);
% plot(0:10, PLSmsep(2,:), 'b-o');
% xlabel('Number of components');
% ylabel('Estimated Mean Squared Prediction Error');
% legend({'PLSR'},'location','NE');

figure
plot(1:10, cumsum(100 * pctVar(2,:)), '-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');

% Weight analysis 
figure
[Xl, Yl, Xs, Ys, beta, pctVar, mse, stats] = plsregress(X, y, nm_pls);
plot(1:size(X, 2), stats.W, '-');
xlabel('Variable');
ylabel('PLS Weight');
legend({'1st Component' '2nd Component' '3rd Component'}, 'location', 'NW');

Y_prd = Xs * Yl'; Y_prd = Y_prd + mean(y);
% true_data = y; predict_data = Y_prd;
true_data = myTbl_matrix(:, end); 
predict_data = Y_prd * std(myTbl_matrix(:, end)) + mean(myTbl_matrix(:, end));
error_pls = error_all(true_data, predict_data);

figure
scatter(true_data, predict_data, 'filled')
tmp_max = max(max([true_data, predict_data])) * 1.2;
tmp_min = min(min([true_data, predict_data])) * 1.2;
axis equal; refline(1,0);
xlim([tmp_min, tmp_max]); ylim([tmp_min, tmp_max]);

[~, large_pls_idx] = maxk(abs(beta(2:end)), nm_comp); % beta(1) is an intercept (not descriptor)--> need to remove
beta2 = beta(2:end);
pls_coeff = beta2(large_pls_idx);
lst_large_pls = pls_coeff';
lst_large_pls_tbl = array2table(lst_large_pls);
lst_large_pls_name = scl_myTbl.Properties.VariableNames(large_pls_idx);
lst_large_pls_tbl.Properties.VariableNames = lst_large_pls_name;
%% Lasso selection of descriptors
k_fold = ipt_prm.k_fold;
y = scl_myTbl_matrix(:,end);
X = scl_myTbl_matrix(:,1:end-1);
[B, FitInfo] = lasso(X,y, 'CV', k_fold);
lassoPlot(B, FitInfo,'PlotType','CV');
legend('show') % Show legend

%Fitinfo includes the index of MinMSE and Index1SE
idx_1se = FitInfo.Index1SE;
B_slct = B(:, idx_1se);
% true_data = y; predict_data = X * B_slct;
% true_data = myTbl_matrix(:, end); 
predict_data = X * B_slct * std(myTbl_matrix(:, end)) + mean(myTbl_matrix(:, end));
error_lasso = error_all(true_data, predict_data);
figure
scatter(true_data, predict_data, 'filled')
tmp_max = max(max([true_data, predict_data])) * 1.2;
tmp_min = min(min([true_data, predict_data])) * 1.2;
axis equal; refline(1,0);
xlim([tmp_min, tmp_max]); ylim([tmp_min, tmp_max]);

large_B_idx = find(B_slct~=0);
lst_large_B = [large_B_idx'; B_slct(large_B_idx)'];
lst_large_B_tbl = array2table(lst_large_B);
lst_large_B_name = scl_myTbl.Properties.VariableNames(large_B_idx);
lst_large_B_tbl.Properties.VariableNames = lst_large_B_name;

% B_slct represents the selected features and coefficients

%% Selection of descriptors by step-wise regression
    tbl = scl_myTbl;
    Criterion = ML_para.stepwise.Criterion;
    PEnter = ML_para.stepwise.PEnter;
%     mdl_stepwise = stepwiselm(tbl,'purequadratic','ResponseVar','PEC');
    mdl_stepwise = stepwiselm(tbl,'constant','Upper','linear','ResponseVar','PEC','Criterion',Criterion,'PEnter',PEnter);
%     mdl_stepwise = stepwiselm(tbl,'linear','ResponseVar','PEC');

    % pureqradratic includs linear and square of each term, interaction
    % includes all product combination
    
    slct_coeff = mdl_stepwise.Coefficients;
    
    true_data = myTbl_matrix(:,end); predict_data = mdl_stepwise.Fitted * std(myTbl_matrix(:, end)) + mean(myTbl_matrix(:, end));;
    error_stpws = error_all(true_data, predict_data);
    figure
    scatter(true_data, predict_data, 'filled')
    tmp_max = max(max([true_data, predict_data])) * 1.2;
tmp_min = min(min([true_data, predict_data])) * 1.2;
axis equal; h = refline(1,0);
xlim([tmp_min, tmp_max]); ylim([tmp_min, tmp_max]);
xlabel("True current [mA/cm^2]")
ylabel("Predicted current [mA/cm^2]")
adjfig
h.Color = 'black';
r_sq_stpws = strcat("R^2 = ",num2str(round(mdl_stepwise.Rsquared.Ordinary,3))); 
annot_stpws = annotation('textbox',[0.3 0.7 0.2 0.15],'FitBoxToText','on','String',r_sq_stpws,'EdgeColor','none');
set(annot_stpws,'FontSize',16,'FontName','Arial')

 
    %%
    tmp1 = mdl_stepwise.VariableNames;
    tmp2 = mdl_stepwise.PredictorNames;
    slct_idx_stpws = [];
    for i=1:size(tmp2,1)
        for j=1:size(tmp1,1)
            if matches(tmp1{j}, tmp2{i})
               slct_idx_stpws = [slct_idx_stpws, j];
            end
        end
    end
    %%
    lst_large_stpws = [slct_idx_stpws; slct_coeff{2:end,1}'];  % Check every time
    lst_large_stpws_tbl = array2table(lst_large_stpws);
    lst_large_stpws_name = tmp2';
    lst_large_stpws_tbl.Properties.VariableNames = lst_large_stpws_name;
    
%% cos similarity mapping and exstraction of descriptors
    dat_matrix = scl_myTbl_matrix;
    scl_myTbl_size = size(scl_myTbl_matrix);
    for i = 1:scl_myTbl_size(2)
        for j = 1:scl_myTbl_size(2)
            cos_sim(i,j) = (dat_matrix(:,i)' * dat_matrix(:,j)) / (norm(dat_matrix(:,i)) * norm(dat_matrix(:,j)));     
        end       
    end
%     cos_sim = dat_matrix' * dat_matrix;
    figure
    h = heatmap(cos_sim);
    h.Colormap = jet;
    
    %% selection of higherst correlation parameters for PEC
    num_large_corr = 11;
    cos_sim_target = cos_sim(end,:);
    [~, large_corr_idx] = maxk(abs(cos_sim_target), num_large_corr);
    lst_large_corr = [large_corr_idx; cos_sim_target(large_corr_idx)];
    lst_large_corr_tbl = array2table(lst_large_corr);
    lst_large_corr_name = scl_myTbl.Properties.VariableNames(large_corr_idx);
    lst_large_corr_tbl.Properties.VariableNames = lst_large_corr_name;
    
    
    %% select the descriptros based on stepwise
        tmp2 = mdl_stepwise.PredictorNames;
         lst_large_stpws_name = tmp2';
    slct_tbl_matrix = ones(myTbl_size(1),length(lst_large_stpws_name));
    for i = 1:length(lst_large_stpws_name)
        tmp_dcp_name = lst_large_stpws_name{i};
        tbl_tmp = scl_myTbl.(tmp_dcp_name);
        slct_tbl_matrix(:,i) = tbl_tmp ;
    end
    
slct_tbl_matrix = [slct_tbl_matrix,scl_myTbl.('PEC')];
    %% fitnlm: nonlinear function fitting (including exp), but specific functional
%  form is necessary.
%% GPR prediction by using all data
    %% training, test index generation
    clear model model_k_fold true_data_trn predict_data_trn 
    num_dat = myTbl_size(1);
%     k_fold = ipt_prm.k_fold;
%     [trn_idx, tst_idx] = idx_generator(num_dat, k_fold);
%     num_dataset = ipt_prm.num_dataset;

    ML_data_trn = slct_tbl_matrix;
%     ML_data_tst = scl_myTbl_matrix(tst_idx(:,num_dataset),:);

%     %% training data and prediction mapping
%     [val_predictions] = GPR_matrix_simple(ML_data_trn);
%     true_data = ML_data_trn(:, end); predict_data_trn = val_predictions;
%     error_trn = error_all(true_data, predict_data_trn);
%     figure
%     scatter(true_data, predict_data_trn, 'filled')
%     tmp_max = max(max([true_data, predict_data_trn])) * 1.2;
%     tmp_min = min(min([true_data, predict_data_trn])) * 1.2;
%     axis equal; refline(1,0);
%     xlim([tmp_min, tmp_max]); ylim([tmp_min, tmp_max]);
    
%% GPR prediction (normal: train and test split)
%% training, test index generation
num_dat = myTbl_size(1);
k_fold = ipt_prm.k_fold;
[trn_idx, tst_idx] = idx_generator(num_dat, k_fold);
%     num_dataset = ipt_prm.num_dataset;
figure
hold on
for j = 1:k_fold
    
    ML_data_trn = slct_tbl_matrix(trn_idx(:,j),:);
    ML_data_tst = slct_tbl_matrix(tst_idx(:,j),:);
    
    %% training data and prediction mapping
    [model, val_predictions, val_RMSE] = GPR_matrix(ML_data_trn);
%     true_data_trn(:,j) = ML_data_trn(:, end); predict_data_trn(:,j) = val_predictions;
    true_data_trn(:,j) = myTbl_matrix(trn_idx(:,j), end); 
    predict_data_trn(:,j) = val_predictions * std(myTbl_matrix(:, end)) + mean(myTbl_matrix(:, end));
    error_trn(j) = error_all(true_data_trn(:,j), predict_data_trn(:,j));
    
    %     figure
    scatter(true_data_trn(:,j), predict_data_trn(:,j),100, 'filled')
    model_k_fold{j} = model;
end
tmp_max = max(max([true_data_trn, predict_data_trn])) * 1.2;
tmp_min = min(min([true_data_trn, predict_data_trn])) * 1.2;
axis equal; h = refline(1,0);
xlim([tmp_min, tmp_max]); ylim([tmp_min, tmp_max]);
xlabel("True current [mA/cm^2]")
ylabel("Prediccted current [mA/cm^2]")

adjfig
hold off
h.Color = 'black';
for i = 1:k_fold
r_sq_train(i) = error_trn(i).R_sq;
end
r_sq_train_ave = mean(r_sq_train);
r_sq_train_str = strcat("R^2 = ",num2str(round(r_sq_train_ave,3))); 
annot = annotation('textbox',[0.55 0.15 0.2 0.15],'FitBoxToText','on','String',r_sq_train_str,'EdgeColor','none');
set(annot,'FontSize',16,'FontName','Arial')
legend({'1st' '2nd' '3rd' '4th' '5th'}, 'location', 'NW');
legend('boxoff')


figure
hold on

for j = 1:k_fold
        ML_data_trn = slct_tbl_matrix(trn_idx(:,j),:);
    ML_data_tst = slct_tbl_matrix(tst_idx(:,j),:);
    model_tmp = model_k_fold{j};
    %% test data and prediction mapping
    tst_predictions = model_tmp.predictFcn(ML_data_tst(:,1:end-1));
%     true_data_tst(:,j) = ML_data_tst(:, end); predict_data_tst(:,j) = tst_predictions;
true_data_tst(:,j) = myTbl_matrix(tst_idx(:,j), end); 
    predict_data_tst(:,j) = tst_predictions * std(myTbl_matrix(:, end)) + mean(myTbl_matrix(:, end));
    error_tst(j) = error_all(true_data_tst(:,j), predict_data_tst(:,j));
    
    scatter(true_data_tst(:,j), predict_data_tst(:,j),100, 'filled')
    
    
    
end
tmp_max = max(max([true_data_tst, predict_data_tst])) * 1.2;
tmp_min = min(min([true_data_tst, predict_data_tst])) * 1.2;
axis equal;revfig; h =refline(1,0);
xlim([tmp_min, tmp_max]); ylim([tmp_min, tmp_max]);
xlabel("True current [mA/cm^2]")
ylabel("Predicted current  [mA/cm^2]")
hold off
% revfig
adjfig
h.Color ='black';
for i = 1:k_fold
r_sq_tst(i) = error_tst(i).R_sq;
end

r_sq_tst_ave = mean(r_sq_tst);
r_sq_tst_str = strcat("R^2 = ",num2str(round(r_sq_tst_ave,3))); 
annot = annotation('textbox',[0.55 0.15 0.2 0.15],'FitBoxToText','on','String',r_sq_tst_str,'EdgeColor','none');
set(annot,'FontSize',16,'FontName','Arial')
legend({'1st' '2nd' '3rd' '4th' '5th' }, 'location', 'NW');
legend('boxoff')



%% Function for summary of the result

function pca_summary

% Assign the directry for ppt file and template
import mlreportgen.ppt.*

pptFile = summary_para.ppt.general.ppt_file ;
pptTemplate = summary_para.ppt.general.ppt_template;
ppt = Presentation(pptFile,pptTemplate);

% Style selection of slide master
masters = getMasterNames(ppt);
masters_str = string(masters);
masters_style = summary_para.ppt.general.masters_style;
masters_idx = matches(masters_str,masters_style);
layoutnames = getLayoutNames(ppt,masters_str(masters_idx));

layoutnames_str = string(layoutnames);
% Add slide for figure(1) and figure(2)

for i = 1:length(TotSlideNum)
    
%     Slides = add(ppt,

end

end




function ppt_generation(ppt_info,varargin)

slide_lists = ppt_info.Properties.RowNames;
Target_calc = slide2show;

for i = 1:size(ppt_info,1)
    target_slide_idx(i,1) = contains(slide_lists{i},Target_calc);
end

target_slide_info = ppt_info(target_slide_idx,:);
ppt_files = target_slide_info.ppt_file;
ppt_template = target_slide_info.ppt_template;

import mlreportgen.ppt.*

for j =1: size(target_slide_info,1)
    % initialization (ppt, master, layout)
    ppts{j,1} = Presentation(ppt_files(j),ppt_template(j));
    %     open(ppts{j,1})
    masters{j,1} = getMasterNames(ppts{j,1});
    masters_str{j,1} = string(masters{j,1});
    masters_target = target_slide_info.Master;
    
    for k = 1:size(masters_str{j,1},2)
        masters_idx(j,k) = matches(masters_str{j,1}(k),masters_target);
    end
    
    layoutnames{j,1} = getLayoutNames(ppts{j,1},masters_str{j,1}(masters_idx(j,:)));
    layoutnames_str{j,1} = string(layoutnames{j,1});
    layoutnames_target = target_slide_info.Layout;
    
    for k = 1:size(layoutnames_str{j,1},2)
        layoutnames_idx(j,k) = matches(layoutnames_str{j,1}(k),layoutnames_target);
    end
    
    slide_format{j,1} = layoutnames_str{j,1}(layoutnames_idx(j,:));
    
    % add the contents
    
    contents_slides{j,1} = add(ppts{j,1},slide_format{j,1});
    
    
    % add title
    SlideTitleObj{j,1} = find(contents_slides{j,1},"SlideTitle");
    SlideTitleRaw(j) = target_slide_info.SlideTitle(j);
    
    replace(SlideTitleObj{j,1}(1),target_slide_info.SlideTitle(j));
    
    
    % add description
    DescriptionObj{j,1} = find(contents_slides{j,1},"Description");
    replace(DescriptionObj{j,1}(1),target_slide_info.Description(j));
    
    %     close(ppts{j,1})
end


   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   







end










































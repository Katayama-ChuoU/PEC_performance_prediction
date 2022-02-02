
load Feature_raman.mat
load Feature_uv_vis.mat
load Feature_xrd.mat
load PEC_target.mat

FeatureData_XRD_nom = FeatureData_XRD;
FeatureData_XRD_nom{:,1:9}=FeatureData_XRD{:,1:9}./FeatureData_XRD{:,5};
FeatureData_Raman_nom = FeatureData_Raman;
FeatureData_Raman_nom{:,1:5}=FeatureData_Raman{:,1:5}./FeatureData_Raman{:,1};


features_tbl = [FeatureData_UV_Vis,FeatureData_Raman_nom,FeatureData_XRD_nom,Target];
% features_tbl(:,18:25) = [];
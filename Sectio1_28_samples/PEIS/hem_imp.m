function [Features_imp_data,freq,freq_log,z_re,z_im,imp_abs,phase_shift]=hem_imp(z_re_all,z_im_all,freq_all,imp_abs_all,phase_shift_all)
%% Data size adjustment 

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

%% Feature extractions
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
         [min_phase_reg{i},min_index_phase_reg{i}] = min(phase_shift_reg{i}); % minimum of phase shift
         min_freq_shift_reg{i} = freq_reg{i}(min_index_phase_reg{i}); % index of minimum of phase shift
         [max_phase_reg{i},max_index_phase_reg{i}] = max(phase_shift_reg{i});
         max_freq_shift_reg{i} = freq_reg{i}(max_index_phase_reg{i});
         
         
         imp_abs_reg{i} = imp_abs(Lmin<=freq_log(:,i) & freq_log(:,i)<=Lmax,:);
         [min_imp_reg{i},min_index_imp_reg{i}] = min(imp_abs_reg{i}); % minimum of phase shift
         min_freq_imp_reg{i} = freq_reg{i}(min_index_imp_reg{i}); % index of minimum of phase shift
         [max_imp_reg{i},max_index_imp_reg{i}] = max(imp_abs_reg{i});
         max_freq_imp_reg{i} = freq_reg{i}(max_index_imp_reg{i});
%                   
         
         
         for j =1:length(freq_all)
             Area_phase_reg{i}(:,j) = trapz(freq_reg{i}(:,j),phase_shift_reg{i}(:,j));
             Area_imp_reg{i}(:,j) = trapz(freq_reg{i}(:,j),imp_abs_reg{i}(:,j));
         end
         
         
         
     end
     
min_phase_reg = cell2mat(min_phase_reg');
min_freq_shift_reg =cell2mat(min_freq_shift_reg');
max_phase_reg = cell2mat(max_phase_reg');
max_freq_shift_reg = cell2mat(max_freq_shift_reg');
Area_phase_reg=cell2mat(Area_phase_reg');

min_imp_reg = cell2mat(min_imp_reg');
min_freq_imp_reg =cell2mat(min_freq_imp_reg');
max_imp_reg = cell2mat(max_imp_reg');
max_freq_imp_reg = cell2mat(max_freq_imp_reg');
Area_imp_reg=cell2mat(Area_imp_reg');

Features_phase = abs(vertcat(min_phase_reg,min_freq_shift_reg,max_phase_reg,max_freq_shift_reg,Area_phase_reg));
Features_imp = abs(vertcat(min_imp_reg,min_freq_imp_reg,max_imp_reg,max_freq_imp_reg,Area_imp_reg));


Features_imp_data = (vertcat(Features_phase,Features_imp));



end



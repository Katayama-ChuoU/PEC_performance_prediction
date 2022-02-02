% for the 28 data about pec correlated with temparature

figure
hold on
for i = 1:7
plot(All_data_bare_hem.PEC.vol(:,i),All_data_bare_hem.PEC.current_density(:,i),'b')

end


for i = 8:14
plot(All_data_bare_hem.PEC.vol(:,i),All_data_bare_hem.PEC.current_density(:,i),'r')
% hold on

end

for i = 15:21
plot(All_data_bare_hem.PEC.vol(:,i),All_data_bare_hem.PEC.current_density(:,i),'color',[0    0.5000         0])
% hold on
end

for i = 22:28
plot(All_data_bare_hem.PEC.vol(:,i),All_data_bare_hem.PEC.current_density(:,i), 'color',[1.0000         0    1.0000])
% hold on
end

hold off


xlabel('Potential /V vs RHE')
ylabel('Current density [mA/cm^2]')

adjfig

%%
temp = [600*ones(7,1);650*ones(7,1);700*ones(7,1);750*ones(7,1);];
figure

scatter(temp,features_tbl.PEC,'filled')

xlabel('Temperature [degree]')
ylabel('Current density [mA/cm^2]')

adjfig
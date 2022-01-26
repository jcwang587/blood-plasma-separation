clear;
clc;

load data
cell_num = 25;

for i=1:cell_num
    plot(cell_tra{i}(:,1),cell_tra{i}(:,3), 'LineWidth', 2);
    hold on
end
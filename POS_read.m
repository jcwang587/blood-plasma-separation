clear;
clc;

path = '\\wsl.localhost\Ubuntu-20.04\home\jcwang\object_in_fluid\output\sim2022-1-10_15-40\';
cell_num = 25;
cell_tra = cell(1,cell_num);
for i = 0:cell_num-1
    tra = [];
    for j = 0:199
        cell_path = [path, 'cell', int2str(i), '_', int2str(j), '.vtk'];
        content = importdata(cell_path, ' ', 5);
        data = content.data;
        pos = mean(data);
        tra = [tra; pos];
        disp([int2str(i),' ',int2str(j)])
    end
    cell_tra{i+1} = tra;
end


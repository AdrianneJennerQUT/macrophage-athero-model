% load files

A = '../output/output0000000';
A2 = '../output/output000000';
A3 = '../output/output00000'; 
A4 = '../output/output0000'; 
B = '.xml';

total_time = 1073;
delta_t_cell = 6;

for tcount = 1:total_time
        clf
    if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    elseif tcount<1001
        K = [A3 num2str(tcount-1,'%d') B];
    else
        K = [A4 num2str(tcount-1,'%d') B];
    end

MCDS = read_MultiCellDS_xml( K , '../output');


cell_types = MCDS.discrete_cells.metadata.type;

mac_C1(tcount) = length(find(cell_types==1));
mac_C2(tcount) = length(find(cell_types==2));

    if isempty(MCDS.discrete_cells.dead_cells==1)
        deadcells(tcount) = 0;
    else
        deadcells(tcount) = length(MCDS.discrete_cells.dead_cells);
    end
end

time = [1:total_time]*delta_t_cell;

figure
subplot(1,3,1)
hold on 
plot(time, mac_C1)
ylabel('C1 macs')

subplot(1,3,2)
hold on
plot(time, mac_C2)
ylabel('C2 macs')

subplot(1,3,3)
hold on
plot(time, deadcells)
ylabel('Dead cells')


% total amount of cells over time of each class 
% amount of internal lipid
% amount of internal lipid in dead cells
% number of dead cells

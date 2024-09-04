% load files

A = '../output/output0000000';
A2 = '../output/output000000';
A3 = '../output/output00000'; 
A4 = '../output/output0000'; 
B = '.xml';

total_time = 1140;
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

%MCDS = read_MultiCellDS_xml( K , '../output');
MCDS = read_MultiCellDS_xml_dbergman( K , '../output');

cell_types = MCDS.discrete_cells.metadata.type;

mac_C1(tcount) = length(find(cell_types==1));
mac_C2(tcount) = length(find(cell_types==2));

    if isempty(MCDS.discrete_cells.dead_cells==1)
        deadcells(tcount) = 0;
    else
        deadcells(tcount) = length(MCDS.discrete_cells.dead_cells);
    end

internal_lipid(tcount) = sum(MCDS.discrete_cells.custom.internalized_total_substrates);

end

%%

time = [1:total_time]*delta_t_cell;

figure
subplot(2,2,1)
hold on 
plot(time, mac_C1)
ylabel('C1 macs')
set(gca,'FontSize',14)

subplot(2,2,2)
hold on
plot(time, mac_C2)
ylabel('C2 macs')
set(gca,'FontSize',14)

subplot(2,2,3)
hold on
plot(time, deadcells)
ylabel('Dead cells')
set(gca,'FontSize',14)

subplot(2,2,4)
hold on
plot(time, internal_lipid./(mac_C1+mac_C2))
ylabel('Avg. internal lipid')
set(gca,'FontSize',14)


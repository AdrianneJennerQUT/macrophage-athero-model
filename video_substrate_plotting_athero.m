% load files

A = '../output/output0000000';
A2 = '../output/output000000';
A3 = '../output/output00000'; 
A4 = '../output/output0000'; 
B = '_microenvironment0.mat';

v = VideoWriter('video1.avi')
v.FrameRate = 10;
open(v)
figure('Position', [50 50 900 700])
hold on 

total_time = 1100;
delta_t_cell = 6;

for tcount = 1:1100
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
    M = read_microenvironment( K ); 
    titles{1} = 'lipid';
    plot_microenvironment( M , titles ); 
    
    lipid_total(tcount) = sum(sum(M.data{1}))*20*20;

    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);

time = [1:total_time]*delta_t_cell;

figure
plot(time,lipid_total)
xlabel('Time')
ylabel('Total lipid')
set(gca,'FontSize',18)

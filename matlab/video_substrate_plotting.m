% Creating movie of substrates
addpath('../output')
var=1;
timetotal =225; % set the total time steps of the data
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000'; 
A4 = 'output0000'; 
B = '.xml';
lipid_index = 7;

v = VideoWriter('video1.avi')
v.FrameRate = 10;
open(v)
figure('Position', [50 50 900 700])
hold on 

for tcount = 2:timetotal
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
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(lipid_index).data(:,:,k) , 20 ) ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
     title( sprintf('%s (%s) at t = %3.2f %s', MCDS.continuum_variables(lipid_index).name , ...
     MCDS.continuum_variables(lipid_index).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units ) ); 
 
 	 lipid(tcount) = sum(sum(MCDS.continuum_variables(lipid_index).data(:,:,k)));
 	 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end
close(v)

%%

figure
hold on
subplot(3,2,1)
plot(virus,'LineWidth',2)
ylabel('Virion')
xlabel('Time (hours)')
set(gca,'FontSize',14)
subplot(3,2,2)
plot(IFN,'LineWidth',2)
xlabel('Time (hours)')
set(gca,'FontSize',14)
ylabel('IFN')
subplot(3,2,3)
plot(chemokine,'LineWidth',2)
xlabel('Time (hours)')
set(gca,'FontSize',14)
ylabel('chemokine')
subplot(3,2,4)
plot(pro_cyto,'LineWidth',2)
xlabel('Time (hours)')
set(gca,'FontSize',14)
ylabel('cytokine')
subplot(3,2,5)
plot(debris,'LineWidth',2)
xlabel('Time (hours)')
set(gca,'FontSize',14)
ylabel('debris')
subplot(3,2,6)
plot(ROS,'LineWidth',2)
xlabel('Time (hours)')
ylabel('ROS')
set(gca,'FontSize',14)
 saveas(gcf,'FZERO.png')

STOP
%%

v = VideoWriter('sub91.avi')
v.FrameRate = 10;
open(v)
figure
hold on 

for tcount = 1:timetotal
    clf
    if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(2).data(:,:,k) , 20 ) ;
    axis image;
    colorbar;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end

close(v)

v = VideoWriter('sub92.avi')
v.FrameRate = 10;
open(v)
figure
hold on 

for tcount = 1:timetotal
     clf
    if tcount<10
        K = [A num2str(tcount,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
   else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(5).data(:,:,k) , 20 ) ;
    axis image
    colorbar 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(5).name , ...
     MCDS.continuum_variables(5).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
    frame = getframe(gcf)
    writeVideo(v,frame);
   

end
close(v)

v = VideoWriter('sub93.avi')
v.FrameRate = 10;
open(v)
figure
hold on 

for tcount = 1:timetotal
    clf
    
    if tcount<10
        K = [A num2str(tcount,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
   else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(7).data(:,:,k) , 20 ) ;
    axis image;
    colorbar ;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(7).name , ...
     MCDS.continuum_variables(7).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end


v = VideoWriter('sub94.avi')
v.FrameRate = 10;
open(v)
figure
hold on 

for tcount = 1:timetotal
    clf
    
    if tcount<10
        K = [A num2str(tcount,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
   else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(9).data(:,:,k) , 20 ) ;
    axis image;
    colorbar ;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(9).name , ...
     MCDS.continuum_variables(9).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end



%%
v = VideoWriter('substrate3and4.avi')
v.FrameRate = 1;
open(v)
figure
hold on 
colormap jet

for tcount = 1:750
    clf
    
    if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(3).data(:,:,k)+MCDS.continuum_variables(4).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    sum(sum((MCDS.continuum_variables(3).data(:,:,k)+MCDS.continuum_variables(4).data(:,:,k)).*20*20))
    axis image;
    colorbar ;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(4).name , ...
     MCDS.continuum_variables(4).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end


%%
%%
v = VideoWriter('sub4.avi')
v.FrameRate = 10;
open(v)
figure
hold on 
colormap jet

for tcount = 1:205
    clf
    
    if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(7).data(:,:,k)+MCDS.continuum_variables(8).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    sum(sum((MCDS.continuum_variables(7).data(:,:,k)+MCDS.continuum_variables(8).data(:,:,k)).*20*20))
    axis image;
    colorbar ;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(4).name , ...
     MCDS.continuum_variables(7).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
    frame = getframe(gcf);
    writeVideo(v,frame);
    

end




%%

% Creating movie of substrates
addpath('../output')

timetotal = 1;
A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

colormap jet
figure
hold on 
tcount =1;
     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(2).data(:,:,k)+MCDS.continuum_variables(2).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 

figure
hold on 
tcount = 30;
colormap jet

     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(2).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
figure
hold on 
tcount = 90;
colormap jet

      if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(2).data(:,:,k)+MCDS.continuum_variables(2).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 

figure
hold on 
tcount = 60;
colormap jet
      if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(3).data(:,:,k)+MCDS.continuum_variables(4).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 
 
 colormap jet
figure
hold on 
tcount =100;
      if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(3).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 

figure
hold on 
tcount = 300;
colormap jet

      if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(3).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
figure
hold on 
tcount = 500;
colormap jet

     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(3).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 

figure
hold on 
tcount = 700;
colormap jet
     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(3).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 

colormap jet
figure
hold on 
tcount =100;
     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(4).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 

figure
hold on 
tcount = 300;
colormap jet

     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(4).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
figure
hold on 
tcount = 500;
colormap jet

     if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(4).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 
 

figure
hold on 
tcount = 700;
colormap jet
      if tcount<10
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<100 
        K = [A2 num2str(tcount,'%d') B]; 
    else 
        K = [A3 num2str(tcount,'%d') B]; 
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    [MCDS.continuum_variables(4).data(:,:,k)] , 150 ,'EdgeColor','none') ;
    axis image;
    colorbar; 
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
     MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units, ... 
     MCDS.mesh.Z_coordinates(k), ...
     MCDS.metadata.spatial_units ) ); 

 
 
 %%
 
% tcount = 2;
 
% tcount = 25;
 

 tcount = 49;
 
 figure
    if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(2).data(:,:,k) , 20 ) ;
    axis image;
    colorbar;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
    title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(2).name , ...%MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units) );%, ... 
     %MCDS.mesh.Z_coordinates(k), ...
    % MCDS.metadata.spatial_units 
     set(gca,'FontSize',13)
    %frame = getframe(gcf);
    %writeVideo(v,frame);
    %%
    
% tcount = 2;
 
% tcount = 25;
 

 tcount = 49;
 
 figure
    if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
    MCDS = read_MultiCellDS_xml(K);
    k = find( MCDS.mesh.Z_coordinates == 0 ); 
    contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), ...
    MCDS.continuum_variables(3).data(:,:,k) , 20 ) ;
    axis image;
    colorbar;
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) ); 
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) ); 
    
  title( sprintf('%s at t = %3.2f %s', MCDS.continuum_variables(3).name , ...%MCDS.continuum_variables(2).units , ...
     MCDS.metadata.current_time , ...
     MCDS.metadata.time_units) );%, ... 
     %MCDS.mesh.Z_coordinates(k), ...
    % MCDS.metadata.spatial_units 
     set(gca,'FontSize',13)
    %frame = getframe(gcf);
    %writeVideo(v,frame);
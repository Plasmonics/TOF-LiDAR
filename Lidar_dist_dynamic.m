%% ------------------------------------
%  TOF LiDAR system signal processing for dynamic platform run
%  Last update: 11/16/2020
%  Author:Hosna Sultana
% For TOF LiDAR timestamp data (time and number of count) for getting range from the maximum probability
% density from the histogram and multiple target range from most prominent peaks

%% ------------------------------------

%  display_by = 1 % plot histogram
%  display_by = 2 % plot min del_t
% 
% if display_by == 1
% elseif display_by == 2  
% end
    
Timestamp_  = xlsread('140a.xlsx');
label = 140;

del_t0 = Timestamp_(:,1);
%time_= Timestamp_(:,1);

del_t_max_W = max(del_t0);       
del_t_min_W = min(del_t0);  

divider = 100;

%divider = 100;
ind = round(length(del_t0)/divider);
A= del_t0';
B = reshape([A,nan(1,ind-mod(numel(A),ind))],[],ind);
B(isinf(B)) = nan;
%fillmissing(x, 'linear'); 
%x = fillmissing(x, 'linear');

del_t = B;

[numRows,numCols] = size(del_t);
%maxcount_ = size(del_t(2));

for i = 1:numCols
    j =  del_t(:,i);
%    
%    del_t(j,i)= del_t(:,i)

nbins= 600000;
%nbins = numRows+2;
   
del_t_max(i) = max(del_t(:,i));       
del_t_min(i) = min(del_t(:,i));   


                        %fixed bin size
                     
%bin_width = 0.001;  %1 ns
%nbins = round((del_t_max- del_t_min)/bin_width);
%[N(i),edges(i)] = histcounts(del_t(:,i)); 


h(i) = histogram(j,nbins,'facecolor',[0.3 0.7 0.2],'edgecolor',[0.3 0.7 0.2]);
%h(i) = histogram(del_t(:,i),nbins,'facecolor',[0.3 0.7 0.2],'edgecolor',[0.3 0.7 0.2]);
[maxcount(i), whichbin(i)] = max(h(i).Values(i));
%title("Maximum of " + maxcount + " occurs in bin " + whichbin + ...
%      " between " + h(i).BinEdges(whichbin) + " and " + h(i).BinEdges(whichbin+1));
 % set(handaxes1, 'Box', 'off')
 grid on
 
 Lt_T_micr_s = 299792458/ 10^6;     % light speed/micro second
%min_dis = Lt_T_micr_s * del_t_min ;  
%max_dis = Lt_T_micr_s * del_t_max ;  



%                for TDC correction

%ave_del_t_min = sum((del_t_min))/numel(del_t_min); 
%Lase_start = ave_del_t_min;
%Laser_Start = 0.1624; % max probability count of laser internal reflection
Laser_Start = 0.150;  %random value for calibration
%Laser_Start = del_t_min;

min_t(i) = (h(i).BinEdges(whichbin(i)) - Laser_Start)/2 ;
max_t(i) = (h(i).BinEdges(whichbin(i)+1) - Laser_Start)/2 ;

%min_t(i) = (h(i).BinEdges(whichbin(i)) - 0.1624)/2 ;
%max_t(i) = (h(i).BinEdges(whichbin(i)+1) - 0.1624)/2 ;


%min_dis(i) = Lt_T_micr_s .* h(i).BinEdges(whichbin(i)) ;  
%max_dis(i) = Lt_T_micr_s .* h(i).BinEdges(whichbin(i)+1) ;
min_dis(i) = Lt_T_micr_s .* min_t(i) ;  
max_dis(i) = Lt_T_micr_s .* max_t(i) ;

%maxcount_ = maxcount_ + maxcount(i);
%i=i+1;
end
%D1 = num2str(h.BinEdges(whichbin),'%100.4d\n');
%D2 = num2str(h.BinEdges(whichbin+1),'%100.4d\n');

%D11(i) = num2str(min_dis(i),'%100.4d\n');
%D22(i) = num2str(max_dis(i),'%100.4d\n');


%time = (1:1:numCols);
Time = 0:0.1:(ind/10)-0.1;
%Time = 0:0.001:(ind/1000)-0.001;
%ind_time = round(length(time_)/divider);
%A_t= time_';
%B_t = reshape([A_t,nan(1,ind_time-mod(numel(A_t),ind_time))],[],ind_time);
%Time = median(B_t);

fontsize = 16;
linewidth = 1;
markersize = 1;

plot (Time', min_dis','.-','color',[1.0 0.56 0.14],'LineWidth',linewidth,'MarkerSize',markersize);
hold on
plot (Time', max_dis','.-','color',[0.0 0.4 0.0],'LineWidth',linewidth,'MarkerSize',markersize);
legend('calculated min distance','calculated max distance','Location','south');
hold off
grid on
annotation('textbox',[0.23 0.88 0.97 0.04],'String',{['Histogram of Lidar data ' num2str(label) ' with nbins = ' num2str(nbins)]},'FitBoxToText','on');

    set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
   % ylabel('Distance in m ')
   % xlabel('time in s')
    %xlabel('time/'num2str(divider)']},'(s)');  
%xlabel( ' Time (s)/ 'Starcat'(num2str(divider));
%xlabel('String',{[' Time (s)/  ' num2str(divider)]});




%                       % del_t max/min plot
% subplot(2,1,1)
% plot(Time',del_t_min','.-','color',[0.2 0.1 0.7],'LineWidth',linewidth,'MarkerSize',markersize);
% legend('minimum value of TOF','Location','south');
% ylabel('TOF ')
% subplot(2,1,2)
% plot(Time',del_t_max','.-','color',[0.8 0.5 0.9],'LineWidth',linewidth,'MarkerSize',markersize);
% legend('maximum value of TOF','Location','south');
% ylabel('TOF ')
%     xlabel('time in s')



    
                          % del_t max/min plot
 
%                           
%  Min_distance = (Lt_T_micr_s.* (del_t_min./2));
%  Max_distance = (Lt_T_micr_s.* (del_t_max./2));
% subplot(2,1,1)
% plot(Time',Min_distance','.-','color',[0.2 0.1 0.7],'LineWidth',linewidth,'MarkerSize',markersize);
% legend('minimum distance from minimum value of TOF','Location','south');
% ylabel('min Distance in m ')
%  xlabel('time in s')
% subplot(2,1,2)
% plot(Time',Max_distance','.-','color',[0.8 0.5 0.9],'LineWidth',linewidth,'MarkerSize',markersize);
% legend('maximum distance from maximum value of TOF','Location','south');
% ylabel('max Distance in m ')
%     xlabel('time in s')
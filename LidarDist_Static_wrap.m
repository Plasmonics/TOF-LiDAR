%% ------------------------------------
%  TOF LiDAR system signal processing for static run
%  Last update: 11/16/2020
%  Author:Hosna Sultana
% For TOF LiDAR timestamp data (time and number of count) for getting range from the maximum probability
% density from the histogram and multiple target range from most prominent peaks

%% ------------------------------------
%Lidar_data is "140a" here

Timestamp_  = xlsread('140a.xlsx');
Timestamp  = Timestamp_(:,1);
label = 140;

nbins = 60000;
del_t = Timestamp;    

%x(isinf(x)) = nan;
%fillmissing(x, 'linear'); 
%x = fillmissing(x, 'linear');
%histogram(del_t67);  

del_t_max = max(del_t);       
del_t_min = min(del_t);  

 bin_width = 0.001;  %1 ns
 nbins = round((del_t_max- del_t_min)/bin_width);

%[nbins,edges] = histcounts(del_t); 

 histogram(del_t,nbins,'facecolor',[0.3 0.7 0.2],'edgecolor',[0.3 0.7 0.2]);

%legend('Timestamp 71','Frontsize','frontsize','Location','northeast');
%legend ( 'Starcat',{[' Timestamp = ' num2str(label), '     calculated  Minimum Distance = ' num2str(D1), ' m, ' , '    calculated  Maximum Distance = ' num2str(D2) ' m ']},'Location','northeast');
%annotation('textbox',[0.23 0.90 0.97 0.024],'String',{['Histogram of Lidar data ' num2str(label) ' with nbins = ' num2str(nbins)]},'FitBoxToText','on');
%set(gca,'fontsize',18);
%grid on          


h = histogram(del_t,nbins,'facecolor',[0.3 0.7 0.2],'edgecolor',[0.3 0.7 0.2]);
[maxcount, whichbin] = max(h.Values);
title("Maximum of " + maxcount + " occurs in bin " + whichbin + ...
      " between " + h.BinEdges(whichbin) + " and " + h.BinEdges(whichbin+1));


%                          %  Distance

table = tabulate(del_t);         
table(:,3) = [];                   
table_max = max(table(:,2));   % how many times one exact data value occurs 
table_min = min(table(:,2));

Lt_T_micr_s = 299792458/ 10^6;     % light speed/micro second


%              without TDC correction
%min_dis = (Lt_T_micr_s * del_t_min)/2 ;  
%max_dis = (Lt_T_micr_s * del_t_max)/2 ;  
%min_dis = Lt_T_micr_s * h.BinEdges(whichbin) ;  
%max_dis = Lt_T_micr_s * h.BinEdges(whichbin+1) ; 


D1 = num2str(h.BinEdges(whichbin),'%100.4d\n');
D2 = num2str(h.BinEdges(whichbin+1),'%100.4d\n');

%                for TDC correction
%Laser_Start = 0.1624;
Laser_Start = del_t_min;
min_t = (h.BinEdges(whichbin) - Laser_Start)/2 ;
max_t = (h.BinEdges(whichbin+1) - Laser_Start)/2 ;
min_dis = Lt_T_micr_s .* min_t ;  
max_dis = Lt_T_micr_s .* max_t ;
D11 = num2str(min_dis,'%100.4d\n');
D22 = num2str(max_dis,'%100.4d\n');


                              % Figure
fontsize = 12;
linewidth = 1;
figure
 histogram(del_t,nbins,'facecolor',[0.3 0.7 0.2],'edgecolor',[0.3 0.7 0.2]);
 legend ( 'Starcat',{[' Timestamp = ' num2str(label), ' with Bin width = ' num2str(bin_width),' (calculated  Minimum Distance = ' num2str(D11), ' m, ' , 'calculated  Maximum Distance = ' num2str(D22) ' m )']},'Location','north');
 set(gca,'fontsize',12);
 
 handaxes1 = axes('position',[0.230,0.425,0.38,0.35]);
    %  fontsize = 7;
    % linewidth = 0.5;
%  markersize = 0.5;
h = histogram(del_t,nbins,'facecolor',[0.3 0.7 0.2],'edgecolor',[0.3 0.7 0.2]);
[maxcount, whichbin] = max(h.Values);
title("Maximum of " + maxcount + " occurs in bin " + whichbin + ...
      " between " + h.BinEdges(whichbin) + " and " + h.BinEdges(whichbin+1));
  set(handaxes1, 'Box', 'off')
 grid on
 set(gca,'fontsize',12);
 xlim([0 0.5]);
 
 
%%
 
 Timestamp_  = xlsread('140a.xlsx');
label = 140;

Timestamp_1hist  = Timestamp_(:,1);


           %to remore too far and too close values
  rowsToDelete_high = (Timestamp_1hist  > 0.5); % get rid of higher TOF values
  Timestamp_1hist (rowsToDelete_high) = [];

%to keep only nearest one target
 rowsToDelete_high_1 = (Timestamp_1hist  > 0.2);
 Timestamp_1hist (rowsToDelete_high_1) = [];
rowsToDelete_low_1 = (Timestamp_1hist  < 0.1499); % get rid of lower TOF values
 Timestamp_1hist (rowsToDelete_low_1) = [];

del_t_1hist = Timestamp_1hist ;

 nbins = 6000000;

 figure
h = histogram(del_t_1hist,nbins);
histfit(del_t_1hist)
% Compute the mean and standard deviation of the actual data.
mu=mean(del_t_1hist);
sigma=std(del_t_1hist);
% Put up lines to indicate the mean, and mean +/- one standard deviation.
line([mu, mu], ylim, 'Color', 'r', 'LineWidth', 1); 
line([mu + sigma, mu + sigma], ylim, 'Color', 'r', 'LineWidth', 1); 
line([mu - sigma, mu - sigma], ylim, 'Color', 'r', 'LineWidth', 1); 
% Make the graph look nice with a grid and labels
grid on
xlabel('difference', 'FontSize', 15);
ylabel('count', 'FontSize', 15);
title('standard deviation of prominent peak', 'FontSize', 15);

yl = ylim; % Get limits of y axis so we can find a nice height for the text labels.
message = sprintf('%.3f  ', mu);
text(mu, 0.95 * yl(2), message, 'Color', 'r', 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
message = sprintf('  %.3f', mu+sigma);
text(mu+sigma, 0.9 * yl(2), message, 'Color', 'r', 'FontSize', 15);
message = sprintf('%.3f  ', mu-sigma);
text(mu-sigma, 0.9 * yl(2), message, 'Color', 'r', 'HorizontalAlignment', 'right', 'FontSize', 15);
 
 
 %%
 Timestamp_  = xlsread('140a.xlsx');
label = 140;

Timestamp  = Timestamp_(:,1);
           %to remore too far and too close values
  rowsToDelete_high = (Timestamp > 0.5); % get rid of higher TOF values
  Timestamp(rowsToDelete_high) = [];
% rowsToDelete_low = (Timestamp < 0.1499); % get rid of lower TOF values
% Timestamp(rowsToDelete_low) = [];

del_t = Timestamp;
% changing the nbins to lower value we can get only range resolution from
% few targets
 nbins = 600;  % change the value
[counts1, binCenters1] = hist(del_t, nbins);

%  plot(binCenters1, counts1, 'r-');
%  
%  legend1 = sprintf('mu = %.3f', mean(del_t));
%  legend({legend1});
%  
                % Pick finding 
  mean_peak_hight =  50;             
 [pks,plocs] = findpeaks(counts1,'MinPeakHeight',mean_peak_hight);                 
%[pks,plocs] = findpeaks([counts1; min(counts1)],'MinPeakProminence',500);

fontsize = 16;
linewidth = 1;
markersize = 1;
figure
plot(binCenters1,counts1)
hold on
%plot(binCenters1(plocs,100), pks, '^r')
plot(binCenters1(plocs), pks, '^r')
hold off
grid on
set(gca,'FontSize',fontsize)
    set(gcf,'Color','w')
    ylabel('counts ')
    %xlabel('time in s')
    annotation('textbox',[0.43 0.88 0.97 0.04],'String',{['Histogram of Lidar data ' num2str(label) '.    with nbins = ' num2str(nbins)]},'FitBoxToText','on');
annotation('textbox',[0.43 0.75 0.97 0.04],'String',{['Mean peak hight ' num2str(mean_peak_hight) '.    number of peaks = ' num2str(length(pks))]},'FitBoxToText','on');


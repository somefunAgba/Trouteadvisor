% Tour Trajectory Advisor Problem (TRAP) using SA (simulated annealing).
%
% Save Cost,
% Save Time,
% Save Data,
% Get Satisfaction,
% Get Profit,
% Get Advised.
%
% (c) July 2018. Team Fuzzy. AI Research Group, EEE Department, FUTA
% (c) oasomefun@ieee.org

% Advisor: Dr E.0. Ogunti
% Credits: GoogleMap
% Inspiration: John Hopfield, George Hinton, J.S. Jang
clearvars
clc;
pause(1);close all;
format shortg;

% fullpath = which(mfilename)
[filepath,name,ext] = fileparts(which(mfilename));
oldpath = cd(filepath);


%% LOAD DATABASE
% Load the driving distance and duration matrix, and the distance route data
dist_mat = load('ng_distdrive_matrix.mat');
duration_mat = load('ng_distduration_matrix.mat');
lineroutes_mat = load('ng_routesroadlatlng.mat');

% load latlon data of most key Nigerian Cities
load('table_ng_caps_lonlat.mat');
map_cityname = [tab_cap.ng_addr]';
loc = [tab_cap.xlngs tab_cap.ylats];

%% Initial Plot
% plot the current tour
figure('Name','SA Tour Trajectory Problem');
fig=gcf;
fig.Units='normalized';
fig.OuterPosition=[0 0 1.2 1];
plot(loc(:, 1), loc(:, 2),'k.','Marker', 's',...
    'Markersize', 6,'MarkerEdgeColor',[0.9 0.9 0.9]);
title(['Tour/Travel Route Trajectory Advisor',''], ...
    'Interpreter','Tex','FontSize',12, 'Color',[0.3 0.5 0.9])
hold on

%% INPUTS
% start_city; % start city input
% toother_city; % other cities input
% if isempty(start_city)
%     start_city = "Akure"; % default
% end
% if isempty(toother_city)
%     toother_city = ["Asaba";"Abakaliki";"Benin City";"Ado Ekiti";...
%     "Enugu"; "Gombe";"Owerri";"Jigawa";"Kaduna"];% default
% end
start_city = "Lagos"; % default

toother_city = ["Asaba";"Abakaliki";"Benin City";"Ilorin";"Ado Ekiti";...
    "Maiduguri";"Enugu";"Gombe";"Owerri";"Jigawa";"Kaduna";...
    "Sokoto";"Akure"];% default


s_check = find(map_cityname == start_city);
if isempty(s_check)
    error('%s \nTry again. \nCity must be a major city in Nigeria.',...
        class(start_city) )
end

num_oc = numel(toother_city);
s_check2 = zeros(1,num_oc);
for i = 1:num_oc
    s_check2(i) = find(map_cityname == toother_city(i));
    if isempty(s_check2)
        error('%s \nTry again. \nAll cities must be a major city in Nigeria.',...
            class(start_city) )
    end
end


%'Plan your Route: By Shortest Average Driving Distance in metres (m)'...
%'or Average Driving Duration in seconds(s)'


%% INITIALIZATIONS
% start index is 1
startIdx = 1;
% number of all the cities from the data map
numCity = length(loc);
num = 1:numCity;
% convert number to 1-D index column vector
num = reshape(num,[numCity,1]);


%% Active city locations
% id for each  selected city in the original data map
act_locs_id = [s_check s_check2]';
% location for the selected cities
act_loc = loc(act_locs_id,:);
% number of selected cities
act_num = length(act_locs_id);
actnum = 1:act_num;
% convert number to 1-D index column vector
actnum = reshape(actnum,[act_num,1]);

%% Distance/Duration matrix for selected cities
act_costmat = zeros(act_num,act_num);
ngCostmat = dist_mat.ng_distdrive_matrix;
for i = 1:act_num
    for j = i:act_num
        act_costmat(i,j) = ngCostmat(act_locs_id(i),act_locs_id(j));
        act_costmat(j,i) = act_costmat(i,j);
    end
end
act2_costmat = zeros(act_num,act_num);
ngCostmat2 = duration_mat.ng_distduration_matrix;
for i = 1:act_num
    for j = i:act_num
        act2_costmat(i,j) = ngCostmat2(act_locs_id(i),act_locs_id(j));
        act2_costmat(j,i) = act2_costmat(i,j);
    end
end


%% SIMULATED ANNEALING
%% initial look-up table
% generating a look-up table to obtain a set of values for dE, change in
% cost energy[where cost => distance or duration] or objective function
% Find typical values of dE
count = 32;
all_dE = zeros(count, 1);
all_dE2 = zeros(count, 1);
for i = 1:count
    % generate a random directional path
    % e.g: id: 7 5 1 8 2 4 3 6
    % means moving from selected city 7 to 5 to 1 to...6
    path = randperm(act_num);
    %% Readjust the random generation
    % Make it a must that the start-city, city 1 starts the path.
    
    % find the index of the start city, city 1
    k = find(path==startIdx);
    % swap city 1 at index 3 with the city 7 at index 1
    % then city 1 will be at index 1 and city 7 at index 3
    path(k) = path(1);
    path(1) = startIdx;
    
    %% Objective function
    % i = [path(2:NumCity) path(1)]
    % ii = (path-1)*NumCity
    % sum(distance(ii+i)) -> cost_energy
    cost_E = sum(act_costmat((path-1)*act_num + [path(2:act_num) path(1)]));
    cost_E2 = sum(act2_costmat((path-1)*act_num + [path(2:act_num) path(1)]));
    %% Move Set:
    new_path = move_sets(act_num,path);
    
    all_dE(i) = abs(cost_E-sum(sum(diff(act_loc([new_path new_path(1)],:))'.^2)));
    all_dE2(i) = abs(cost_E2-sum(sum(diff(act_loc([new_path new_path(1)],:))'.^2)));
end

dE = max(all_dE);dE2 = max(all_dE2);
fprintf('1. Initial Energy Cost = %f\n\n',cost_E);

% get city names of the the selected cities from the data map
tour_city = map_cityname(act_locs_id);
% initial tour and back to start point
out = [path path(1)];
%% First Plot: Show the active cities on the map
% annotate selected city latlon locations on the map
plot(act_loc(out(:), 1), act_loc(out(:), 2),'k.','Marker', 's',...
    'Markersize', 6,'MarkerEdgeColor','m');
h = line(act_loc(out(:), 1), act_loc(out(:), 2),'LineWidth',1.5,...
    'LineStyle','-.', 'Color',[0.8 0.8 0.8]);
% annotate selected city visit-numbers and names in the graphic window
tx = text(act_loc(path(:),1)-0.09, act_loc(path(:),2)+0.07,int2str(actnum(:)), ...
    'Interpreter','Tex','FontSize', 9, 'Color',[0.1 0.1 0.1]);
cx = text(loc(:, 1), loc(:, 2)-0.15, map_cityname(:),...
    'Interpreter','tex','FontSize',8,'FontWeight','bold',...
    'Color',[0.85 0.85 0.85]);
delete(cx(act_locs_id));
acx = text(act_loc(:, 1), act_loc(:, 2)-0.15, tour_city(:),...
    'Interpreter','tex','FontSize',8,'FontWeight','bold',...
    'Color',[0.3 0.5 0.9]);
hold off
grid on; grid minor

%% ANNEALING SCHEDULE
% Choose the temperature to be large enough
an.temp = 10*dE; % modified to be very small
an.temp2 = 10*dE2; % modified to be very small
an.MaxTrialN = act_num*100; % Max. # of trials or iterations at a temperature
an.MaxAcceptN = act_num*10; % Max. # of acceptances at a temperature
an.TempRatio = 0.5; % Temperature decrease ratio
an.StopTolerance = an.TempRatio/1e6;	% Stopping tolerance
an.minE = inf;	% Initial value for min. dist_energy
an.maxE = -1;	% Initial value for max. dist_energy
an.dist_energy = cost_E;an.dist_energy2 = cost_E2;
% Major annealing loop
it = 0;
while (an.temp > 1e-8) || (an.maxE - an.minE)/an.maxE > an.StopTolerance
    an.TrialN = 0;	% Number of trial moves
    an.AcceptN = 0; % Number of actual moves
    %% TRIALS
    [an,path] = annealing(an,act_num,act_costmat,path);
    [an2,path] = annealing(an,act_num,act2_costmat,path);
    
    %% INFO.
    % Final Temperature
    fprintf('Temperature = %f\n', an.temp);
    % Final Path- % tmp = sprintf('%d ',path);
    fprintf('Tour Path = %s\n', num2str(path));
    fprintf('Final dist_energy, E = %f\n', an.dist_energy);
    fprintf('[min-E max-E] = [%f %f]\n', an.minE, an.maxE);
    fprintf('[NumberAccepted NumberOfTrials ScheduleIterations]= [%d %d %d]\n\n',...
        an.AcceptN, an.TrialN, it);
    
    % Lower the temperature
    an.temp = an.temp*an.TempRatio;
    
    curr_out = [path path(1)];

    %% SESSION BEST-SOLUTION
    minDist = an.dist_energy;
    minTime = an.dist_energy2;

    best_path = path;
    best_out = [best_path best_path(1)];
    
    %% UPDATE PLOT
    out = best_out;
    
    % Outline Driving Routes
    actv_outlinepath = [];
    for i = 1:act_num
        j = i + 1;
        actv_outlinepath = [[actv_outlinepath]; 
        lineroutes_mat.ng_routesroadlatlng{act_locs_id(out(i)),act_locs_id(out(j))}];
    end
    if exist('outl', 'var')
        delete(outl);
    end
    outl = line(actv_outlinepath(:, 1)', actv_outlinepath(:, 2)',...
        'Color',[0.5 0.9 0.36],'LineStyle', '-.','LineWidth', 1.5);
    %     [0.25 1 0.50][0.6 0.89 0.36]
    
    % Display Cost
%     if cost =='m'
        cost_str_d = ['Optimal Trajectory Cost : ', num2str(an.dist_energy/1000), ' km'];
%     else
        tmp = an.dist_energy2/86400;
        DD = floor(tmp);
        tmp= (tmp - DD)*24;
        HH = floor(tmp);
        tmp = (tmp - HH)*60;
        MM = floor(tmp);
        cost_str_t = ['Optimal Trajectory Cost : ', num2str(DD),' days ',...
            num2str(HH),' hours ', num2str(MM), ' minutes '];
%     end
    
    % Display Cost as Legend
    lgd = legend({cost_str_d, cost_str_t}, 'Location','north','FontSize',...
        9,'FontWeight','bold','Interpreter','Tex','Box','off',...
        'TextColor',[0.9,0.42,0.45]);
    
    % Update Route Map
    set(h, 'xdata', act_loc(out(:), 1), 'ydata', act_loc(out(:), 2));
    delete(tx);
    tx = text(act_loc(best_path(:),1)-0.15, act_loc(best_path(:),2)+0.15,...
        int2str(actnum(:)), ...
        'Interpreter','Tex','FontSize', 8,'FontWeight','bold',...
        'Color',[0.9,0.42,0.45]);
    drawnow limitrate;
end

% Create directory to store views of best route
ulxy = [act_loc(best_path(act_num),1) act_loc(best_path(act_num),2)];
ulid = find(act_loc==[ulxy(1) ulxy(2)]);
if (ulid(2)-ulid(1)) == act_num
    ulidx = ulid(1);
    last_city = tour_city(ulidx);
    
    namef = strcat('trp-saves');
    [status, msg, msgID] = mkdir (namef);
    if status == 0
        disp(msg)
    end
    % Save File
    oldPath = cd(namef);
    namefile = strcat('trp-',num2str(start_city),...
        '-',last_city,num2str(act_num));
    fPath = strcat(oldPath,'\',namef,'\',namefile,'.png');
    if exist(fPath,'file')
        namefile = strcat(namefile,'-n');
    end
    saveas(gcf,strcat(namefile,'.png'))
    cd(oldPath)
end



% zoom, pan
% hCM = uicontextmenu;
% hMenu = uimenu('Parent',hCM,'Label','Pan',...
% 'Callback','pan(gcbf,''on'')'); %#ok<NASGU>
% hZoom = zoom(gcf);
% hZoom.UIContextMenu = hCM;
zoom('on')

% hMenu = uimenu('Parent',hCM,'Label','Zoom',...
%     'Callback','zoom(gcbf,''on'')');
% hPan = pan(gcf);
% hPan.UIContextMenu = hCM;
% pan('on')


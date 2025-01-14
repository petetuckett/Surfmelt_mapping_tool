%% POST-PROCESSING SCRIPT FOR GOOGLE EARTH ENGINE GEOJSON FILES

% Use this script to convert GEOJSON output files from GEE into useable shapefiles.
% This script performs various filtering and merging steps, and adds metadata to lake polygons.

% (1) Read GeoJSON file
% (2) Loop through time windows and output shapefiles for each window for each catchment
% (3) Combine catchment specific shapefiles into one big shapefile per time window
% (4) Merges lake shapefiles, cuts out islands, and assigns new metadata (incl elevation and grounded)

% File pathways will need changing by the user. A High Performace Computing (HPC) platform may be required to run the script.

% Pete Tuckett 2024
% pete.tuckett@york.ac.uk, petetuckett47@gmail.com

% See this paper for details of the mapping methodology:
% Tuckett, P. A., Ely, J. C., Sole, A. J., Lea, J. M., Livingstone, S. J., Jones, J. M., & van Wessem, J. M. (2021). 
% Automated mapping of the seasonal evolution of surface meltwater and its links to climate on the Amery Ice Shelf, Antarctica. The Cryosphere. 1-35.


%% User Parameters

addpath(genpath(['/mnt/fastdata/ggp19pat/Lake_Mapping/Post_processing']));

Region = 'AIS_Q2'; % NOTE - region folders need to exist prior to running script
Satellite = 'LST';
% window_length_days = 15; % number of days between images as set in GEE -----------------> Need to change this method if using bimonthly method rather than set number of window days
geojson_filename = 'LST_AIS_Q2_20150701_20160630.geojson';
start_date = '20150701'; % Input start date (same as GEE)
end_date = '20160630'; % Input end date (same as GEE)
Min_lake_area = 0.0018; % Minimum lake area size for area filter (0.0002 for S2)

Processing_stage = 4;
% 1 = First stage only
% 2 = Stages 1 and 2
% 3 = Stages 1, 2 and 3
% 4 = All Stages (1-4)


%% (1) Read GeoJSON file

disp('STAGE 1: Reading GeoJSON file')

% Set Input folder that contains the GeoJSON file
%inputPath = 'C:\Users\Peter\Documents\PhD\Earth_Engine\geojson_files\';
inputPath = (['../geojson_files/Continent_wide_results/' Region '/']);

% Set the Output folder in which shapefiles will be saved
%outputPath = 'C:\Users\Peter\Documents\PhD\Earth_Engine\Shapefile_outputs\Rennick\Ren_Jan2018\';
mkdir (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite]);
mkdir (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite '/Individial_tiles']);
outputPath = (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite '/Individial_tiles/']);

disp('Loading geoJSON file');

% Finds filepath for geojson file
filename = strcat(inputPath,geojson_filename);

input_file = fileread(filename);  % Reads geoJSON file
input = jsondecode(input_file);   % Decodes JSON file
disp('geoJSON loaded :-)')

% Converts start and end dates to datenum format
start_datenum = datenum(start_date,'yyyymmdd');
end_datenum = datenum(end_date,'yyyymmdd');

%% (2) Run a loop through all the windows to extract shapefiles

if Processing_stage >= 2

disp('STAGE 2: Looping through time windows to extract raw shapefiles')
warning('off')

numWindows = size(input.features);

% Loops through all features to calculate a Window Visibility Score based on PercVis and PercVal statistics, adds as property to each feature
for p = 1:numWindows(1,1)
    PercVis = input.features(p).properties.PercVis;
    PercVal = input.features(p).properties.Perc_values;
    FNVis = fieldnames(PercVis);
    FNVal = fieldnames(PercVal);
    NumCommonfields = ismember(FNVal, FNVis);
    PercVisMatchedFields = struct();
        % Loop to get rid of PercVis values from images that are not included within PercVal
        for K = 1 : length(NumCommonfields)
          if NumCommonfields(K); PercVisMatchedFields.(FNVal{K}) = PercVis.(FNVal{K}); end
        end
    CombinedStructures = [PercVisMatchedFields, PercVal];
    Combtable = struct2table(CombinedStructures);
    Combarray = table2array(Combtable);
    Fractions = Combarray(2,:)/100; % Divide percentages by 100
    Combarray(2,:) = Fractions;
    VisScores = prod(Combarray); % Multiply Vis score by PercVal (e.g. 90% visibiity * 0.3)
    WindowVisScore = sum(VisScores);
    
    % Add the final visibility score as a property to the feature
    input.features(p).properties.CompositeVisScore = WindowVisScore;
   
end

% Create a table of metadata for each tile for each time window - used in stage 4 when re-assigning metadata to merged lake shapefiles
Tile_metadata = [];
for n = 1:numWindows(1,1)
    dateStart = string(input.features(n).properties.dateStart);
    dateEnd = string(input.features(n).properties.dateEnd);
    tile = string(input.features(n).properties.catchment);
    numImages = input.features(n).properties.numImages;
    CompositeVisScore = input.features(n).properties.CompositeVisScore;
    StackedVisScore = input.features(n).properties.PercentVisible;
    TableInfo = table(dateStart, dateEnd, tile, numImages, CompositeVisScore, StackedVisScore);
    % TableInfo = table(dateStart, dateEnd, numImages, CompositeVisScore, StackedVisScore);

    Tile_metadata = cat(1,[Tile_metadata;TableInfo]);
end
% Convert table to a structure
Tile_metadata = table2struct(Tile_metadata);

% Start big loop to go through each feature, will then extract lake geometries
for n = 1:numWindows(1,1)
    output = [];
    a = 0;
    master_datenum = datenum(num2str(input.features(n).properties.dateStart,'yyyy-mm-dd')); % In datenum format
    slave_datenum = datenum(num2str(input.features(n).properties.dateEnd,'yyyy-mm-dd'));
    
    % dateEnd = datestr(datenum(slave_datenum),'yyyy-mm-dd'); % Creates string version of dateEnd
    
    ds = day(datetime(input.features(n).properties.dateStart,'InputFormat','yyyy-MM-dd')); % ds (day start)
    ms = month(datetime(input.features(n).properties.dateStart,'InputFormat','yyyy-MM-dd')); % ms (month start)
    ys = year(datetime(input.features(n).properties.dateStart,'InputFormat','yyyy-MM-dd')); % ys (year start)
    
    de = day(datetime(slave_datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyyMMdd')); % de (day end)
    me = month(datetime(slave_datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyyMMdd')); % me (month end)
    ye = year(datetime(slave_datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyyMMdd')); % ye (year end)
    
    catchment = input.features(n).properties.catchment; % Comment out if running on a single ROI
    
    % Ensure file names are always in the format yyyymmdd_yyyymmdd
    ds = sprintf('%02d', ds);
    ms = sprintf('%02d', ms);
    de = sprintf('%02d', de);
    me = sprintf('%02d', me);
    
    % Create output filename based on master and slave dates
    %outputFilePath = strcat(outputPath,num2str(ys),num2str(ms),num2str(ds),'_',num2str(ye),num2str(me),num2str(de),'.shp'); % If running on single ROI
    outputFilePath = strcat(outputPath,Satellite,'_',num2str(ys),num2str(ms),num2str(ds),'_',num2str(ye),num2str(me),num2str(de),'_',Region,'__',catchment,'.shp'); % If running across multiple catchments
    
    
    % e.g. outputFilePath = ...Window_shapefiles/20200412_20200428.shp
    
    disp(['Writing lake shapefiles for...  LST_' num2str(ys),num2str(ms),num2str(ds) '_' num2str(ye),num2str(me),num2str(de) '_' Region '__' catchment]);
    % disp(['Writing lake shapefiles for...  L8_' num2str(ys),num2str(ms),num2str(ds) '_' num2str(ye),num2str(me),num2str(de) '_' Region]);

    
   % Extract the coordinates and metadata for each individual lake
    if strcmp(input.features(n).geometry.type,'MultiPolygon')
        coordinates = squeeze(input.features(n).geometry.coordinates);
        if iscell(coordinates)
             for m = 1:length(coordinates)
                if length(coordinates{m}(:,1,1)) == 1
                    a = a+1;
                    [X,Y] = polarstereo_fwd(coordinates{m}(1,:,2), coordinates{m}(1,:,1), 6378137.0, 0.08181919, -71, 0);
                    output(a).X = X;
                    output(a).Y = Y;
                    output(a).Geometry = 'Polygon';
                    output(a).dateStart = input.features(n).properties.dateStart;
                    output(a).dateEnd = input.features(n).properties.dateEnd;
                    output(a).catchment = input.features(n).properties.catchment;
                    output(a).numImages = input.features(n).properties.numImages;
                    output(a).windowDays = input.features(n).properties.windowDays;
                    output(a).CompositeVisScore = input.features(n).properties.CompositeVisScore;
                    output(a).StackedVisScore = input.features(n).properties.PercentVisible;
                    output(a).lakeArea = polyarea(X,Y)*1e-6;
                    [output(a).centroidX,output(a).centroidY] = centroid(polyshape(X,Y));
                else
                    if length(coordinates{m}(1,:,1)) == 1
                        for k = 1:length(coordinates{m}(:,1,1))
                            a = a+1;
                            [X,Y] = polarstereo_fwd(coordinates{m}{k}(:,2), coordinates{m}{k}(:,1), 6378137.0, 0.08181919, -71, 0);
                            output(a).X = X;
                            output(a).Y = Y;
                            output(a).Geometry = 'Polygon';
                            output(a).dateStart = input.features(n).properties.dateStart;     
                            output(a).dateEnd = input.features(n).properties.dateEnd;
                            output(a).catchment = input.features(n).properties.catchment;
                            output(a).numImages = input.features(n).properties.numImages;
                            output(a).windowDays = input.features(n).properties.windowDays;
                            output(a).CompositeVisScore = input.features(n).properties.CompositeVisScore;
                            output(a).StackedVisScore = input.features(n).properties.PercentVisible;
                            output(a).lakeArea = polyarea(X,Y)*1e-6;
                            [output(a).centroidX,output(a).centroidY] = centroid(polyshape(X,Y));
                        end
                    else
                        for k = 1:length(coordinates{m}(:,1,1))
                            a = a+1;
                            [X,Y] = polarstereo_fwd(coordinates{m}(k,:,2), coordinates{m}(k,:,1), 6378137.0, 0.08181919, -71, 0);
                            output(a).X = X;
                            output(a).Y = Y;
                            output(a).Geometry = 'Polygon';
                            output(a).dateStart = input.features(n).properties.dateStart;
                            output(a).dateEnd = input.features(n).properties.dateEnd;
                            output(a).catchment = input.features(n).properties.catchment;
                            output(a).numImages = input.features(n).properties.numImages;
                            output(a).windowDays = input.features(n).properties.windowDays;
                            output(a).CompositeVisScore = input.features(n).properties.CompositeVisScore;
                            output(a).StackedVisScore = input.features(n).properties.PercentVisible;
                            output(a).lakeArea = polyarea(X,Y)*1e-6;
                            [output(a).centroidX,output(a).centroidY] = centroid(polyshape(X,Y));
                        end
                    end
                end
             end
        else
                for k = 1:length(coordinates(:,1,1))% :,1,1
                    a = a+1;
                    [X,Y] = polarstereo_fwd(coordinates(k,:,2), coordinates(k,:,1), 6378137.0, 0.08181919, -71, 0);
                    output(a).X = X;
                    output(a).Y = Y;
                    output(a).Geometry = 'Polygon';
                    output(a).dateStart = input.features(n).properties.dateStart;
                    output(a).dateEnd = input.features(n).properties.dateEnd;
                    output(a).catchment = input.features(n).properties.catchment;
                    output(a).numImages = input.features(n).properties.numImages;
                    output(a).windowDays = input.features(n).properties.windowDays;
                    output(a).CompositeVisScore = input.features(n).properties.CompositeVisScore;
                    output(a).StackedVisScore = input.features(n).properties.PercentVisible;
                    output(a).lakeArea = polyarea(X,Y)*1e-6;
                    [output(a).centroidX,output(a).centroidY] = centroid(polyshape(X,Y));
                end
        end
       
    elseif strcmp(input.features(n).geometry.type,'Polygon')
        coordinates = squeeze(input.features(n).geometry.coordinates);
        if iscell(coordinates)
            if length(coordinates(1,:,1)) == 1
                for k = 1:length(coordinates(:,1,1))
                    a = a+1;
                    [X,Y] = polarstereo_fwd(coordinates{k}(:,2), coordinates{k}(:,1), 6378137.0, 0.08181919, -71, 0);
                    output(a).X = X;
                    output(a).Y = Y;
                    output(a).Geometry = 'Polygon';
                    output(a).dateStart = input.features(n).properties.dateStart;     
                    output(a).dateEnd = input.features(n).properties.dateEnd;
                    output(a).catchment = input.features(n).properties.catchment;
                    output(a).numImages = input.features(n).properties.numImages;
                    output(a).windowDays = input.features(n).properties.windowDays;
                    output(a).CompositeVisScore = input.features(n).properties.CompositeVisScore;
                    output(a).StackedVisScore = input.features(n).properties.PercentVisible;
                    output(a).lakeArea = polyarea(X,Y)*1e-6;
                    [output(a).centroidX,output(a).centroidY] = centroid(polyshape(X,Y));
                end
            end
        else
        a = a+1;
        [X,Y] = polarstereo_fwd(coordinates(:,2), coordinates(:,1), 6378137.0, 0.08181919, -71, 0);
        output(a).X = X;
        output(a).Y = Y;
        output(a).Geometry = 'Polygon';
        output(a).dateStart = input.features(n).properties.dateStart;
        output(a).dateEnd = input.features(n).properties.dateEnd;
        output(a).catchment = input.features(n).properties.catchment;
        output(a).numImages = input.features(n).properties.numImages;
        output(a).windowDays = input.features(n).properties.windowDays;
        output(a).CompositeVisScore = input.features(n).properties.CompositeVisScore;
        output(a).StackedVisScore = input.features(n).properties.PercentVisible;
        output(a).lakeArea = polyarea(X,Y)*1e-6;
        [output(a).centroidX,output(a).centroidY] = centroid(polyshape(X,Y));  
        end
    end
    
    % Write shapefile and save to output folder
    if length(output)>0
        shapewrite(output,outputFilePath);
    end
    
    clearvars output master_datenum slave_datenum
end

else
    disp('Not running Stage 2');
end

%% (3) Combine catchment specific shapefiles into one big shapefile per time window

if Processing_stage >= 3
   
disp('STAGE 3: Combining catchment specific shapefiles')

mkdir (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite '/Combined']);
CombinedOutputPath = (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite '/Combined/']);

Shapefile_folder = outputPath;
ShapeFiles = dir(fullfile(Shapefile_folder,'*.shp'));

% Create an array of datenums which indicate the start date of each time window
%Window_start_dates_old = start_datenum:window_length_days:end_datenum;

% New Method to generate array of datenums to indicate start of time windows...
startdateoutput = [];

for g = 1:numWindows(1,1)
    startDate = datenum(num2str(input.features(g).properties.dateStart,'yyyy-mm-dd'));
    startdateoutput = cat(1,startdateoutput,startDate);
end

Window_start_dates = unique(startdateoutput);

%Loop through each window start date...
  for W = 1:length(Window_start_dates)
    
    Windowstart = Window_start_dates(W);
    WindowstartDatetime =  datestr(datetime(Windowstart,'ConvertFrom','datenum'));

    Allfiles = [];
    disp(['Writing combined lake shapefile for time window starting... ' WindowstartDatetime])
    
    % Loop through each catchment specific shapefile and combine together if the window start date matches the start date specified in the W loop
    for e = 1:length(ShapeFiles)
        basefilename = ShapeFiles(e).name;
        fullfilename = fullfile(Shapefile_folder, basefilename);
        file = shaperead(fullfilename);
        
        if  datenum(num2str(file(1).dateStart),'yyyy-mm-dd') == Windowstart
            filetmp = file;
            Allfiles = cat(1,[Allfiles;filetmp(:)]);
            datesfilename = strcat(basefilename(1:21),'_',Region,'_comb.shp'); % YYYYMMDD_YYYYMMDD_xxxxxx_all
            
            % Set new output pathway and write combined shapefile
            NewOutputPath = strcat(CombinedOutputPath,datesfilename);
            if length(Allfiles)>0
                shapewrite(Allfiles,NewOutputPath);
            end
        end
    end
    
    clear Allfiles
  end
  
else
    disp('Not running Stage 3')
end

%% (4) Merges lake shapefiles, and assigns new metadata (incl elevation and grounded)

if Processing_stage == 4

warning('off')
disp('STAGE 4: Fixing bugs in raw shapefile outputs and writing new processed shapefiles')

Shapefile_folder = CombinedOutputPath;
%Shapefile_folder = 'C:\Users\Peter\Documents\PhD\Earth_Engine\Shapefile_outputs\Results_v3\Amery\Amery_20170101_20171231\LST\Combined'
ShapeFiles = dir(fullfile(Shapefile_folder,'*.shp'));

mkdir (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite '/Filtered_shapefiles']);
FilteredOutputPath = (['../Shapefile_outputs/Continent_wide_results_v2/' Region '/' Region '_' start_date '_' end_date '/' Satellite '/Filtered_shapefiles/']);

tic;

% % Load datasets needed for metadata
% % Load REMA geotiff for elevation extraction
% disp('Loading REMA 100 m DEM...');
% [elev,R] = geotiffread('../Datasets/REMA/REMA_100m_dem.tif');
% 
% % Load REMA slope map geotiff for slope extraction
% disp('Loading Slope 100 m Map...');
% [slope,A] = geotiffread('../Datasets/REMA/REMA_100m_slope.tif');
% 
% % Load Measures velocity geotiff for velocity extraction
% disp('Loading Measures Velocity 450 m...');
% [velocity,Z] = geotiffread('../Datasets/Measures_velocity/Regional_geotiffs/Amery_velocity.tif');
% 
% % Load Euclidean Distance geotiff for extraction of lake distances to grounding line
% disp('Loading Euc_Distance to Grounding Line 200 m...');
% [Euc_distance,B] = geotiffread('../Datasets/Depoorter_grounding_line/GL_Euc_Distance_200m.tif');

% Load Bedrock Distance geotiff for extraction of lake distances to exposed
% bedrock (or coastline) - based off ice mask made in GEE
% disp('Loading Euc_Distance to Bedrock 60 m...');
% [Bedrock_distance,C] = geotiffread('C:\Users\Peter\Documents\PhD\Earth_Engine\IceMasks\Bedrock_Euc_Distance_Amery_60m.tif');

% Load Grounding line dataset (Depoorter: https://doi.pangaea.de/10.1594/PANGAEA.819147)
disp('Loading Grounding line dataset (Depoorter)...');
Groundingline = shaperead('../Datasets/Antarctica_masks/scripps_antarctica_polygons_v1.shp');
% AJS: Identify the grounded ice sheet
[val, idx] = max(cat(1,Groundingline(:).Area_km2));
Groundingline = Groundingline(idx);
toc;

disp(['------------------------------------------------------------']);

%% Loop through shapefiles 
for S = 1:length(ShapeFiles)%3:3(ShapeFiles) Need to run 1:3...
    
    basefilename = ShapeFiles(S).name;
    outputfilename = strcat(basefilename(1:28),'_filt.shp');
    disp(['Starting ' outputfilename '...']);
    fullfilename = fullfile(Shapefile_folder, basefilename);
    file = shaperead(fullfilename);
    OutputFilePath = strcat(FilteredOutputPath,outputfilename);
    
    
    % Process polygon topology
    disp(['- Simplifying and merging polygons...']);
    % - Concatenate all vertices
    Xall = [];
    Yall = [];
    for n = 1:length(file)
        % Apply lake area filter
        if file(n).lakeArea >= Min_lake_area
        Xtmp = file(n).X;
        Ytmp = file(n).Y;
        Xall = cat(1,[Xall;Xtmp(:)]);
        Yall = cat(1,[Yall;Ytmp(:)]);
        end
    end
    % - Split polygons based on occurrence of NaNs
    [Ycells,Xcells] = polysplit(Yall,Xall);
    % - Convert split polygons to polyshapes, retaining nested polygons
    %  ('holes')
    shapes = polybuffer(polyshape(Xcells,Ycells),1);
    % - Simplify and buffer the polyshapes
    new_output = polybuffer(simplify(union(shapes,'KeepCollinearPoints',false)),-1); % Simplify first??
    toc;
    
    % Make new polyshapes of reach region (i.e. each lake)
    lake1 = regions(new_output);
    
    if length(lake1)>0
        
    disp('- Calculating lake elevation and grounding status...');
    for n = 1:length(lake1)
        lake(n).X = lake1(n).Vertices(:,1);
        lake(n).Y = lake1(n).Vertices(:,2);
        lake(n).Geometry = 'Polygon';
        lake(n).dateStart = file(n).dateStart; % Works because all the same... (length of file is different to length of lake)
        lake(n).dateEnd = file(n).dateEnd;
        lake(n).numImages = file(n).numImages;
        lake(n).windowDays = file(n).windowDays;
        lake(n).lakeArea = area(lake1(n))*1e-6;
        [lake(n).centroidX,lake(n).centroidY] = centroid(lake1(n));
        %lake(n).lakeArea = polyarea(lake1(n).Vertices(:,1),lake1(n).Vertices(:,2))*1e-6;
        %[lake(n).centroidX,lake(n).centroidY] = centroid(polyshape(lake1(n).Vertices(:,1),lake1(n).Vertices(:,2)));
        %lake(n).PercVis = file(n)
    end
    
    % Load Antarctic grid shapefile and assign newly merged lakes to a catchment tile based on X and Y Centroids
    Antarctic_grid = shaperead('../Datasets/ROIs/Antarctica_ROIs.shp');
    Antarctic_grid_v2 = struct2table(Antarctic_grid);
    
    %testout = lake(all(~cellfun(@isnan,struct2cell(lake)))); % Ensures lakes with blank tile values are removed prior to running the loop (otherwise n exceeds the array length if removed within the loop)
    
        for n = 1:length(lake)
            cenX = lake(n).centroidX;
            cenY = lake(n).centroidY;
            idx = cellfun(@(X,Y) inpolygon(cenX,cenY,X,Y),Antarctic_grid_v2.X, Antarctic_grid_v2.Y);
            tile = Antarctic_grid_v2.tile(idx);
                if isempty(tile)
                   lake(n).tile = string('Nan'); % Assigns NaN if tile is missing
                else
                    lake(n).tile = string(tile); %string()?
                end
        end
     
        lakeTab = struct2table(lake);
        lakeTab(find(lakeTab.tile == 'Nan'),:) = []; % Removes rows where the tile is NaN
        lake = table2struct(lakeTab); % Converts back to a structure
    
        % Now need to re-assign numImages and PercentIceVis based on new lakes
        for p = 1:length(Tile_metadata) %1:length(Tile_metadata)
            tile_id = Tile_metadata(p).tile;
            dateStart = Tile_metadata(p).dateStart;
                for n = 1:length(lake)
                    if  lake(n).tile == tile_id && lake(n).dateStart == dateStart
                        lake(n).numImages = double(Tile_metadata(p).numImages);
                        lake(n).CompositeVisScore = double(Tile_metadata(p).CompositeVisScore);
                        lake(n).StackedVisScore = double(Tile_metadata(p).StackedVisScore);
                    end
                end
        end
                    
%         % Interpolates REMA elevation for the centre of each lake
%         lakeElev = mapinterp(elev,R,[lake(:).centroidX],[lake(:).centroidY]);
%         for m = 1:length(lakeElev)
%             lake(m).lakeElev = double(lakeElev(m));
%         end
%         
%         % Interpolates REMA slope for the centre of each lake
%         lakeSlope = mapinterp(slope,A,[lake(:).centroidX],[lake(:).centroidY]);
%         for m = 1:length(lakeSlope)
%             lake(m).Slope = double(lakeSlope(m));
%         end
%         
%         % Interpolates Measures Velocity for the centre of each lake
%         lakeVelocity = mapinterp(velocity,Z,[lake(:).centroidX],[lake(:).centroidY]);
%         for m = 1:length(lakeVelocity)
%             lake(m).Velocity = double(lakeVelocity(m));
%         end
%         
%         % Interpolates Distance to grounding line for the centre of each lake
%         lakeGLDist = mapinterp(Euc_distance,B,[lake(:).centroidX],[lake(:).centroidY]);
%         for m = 1:length(lakeGLDist)
%             lake(m).GL_Distance = double(lakeGLDist(m));
%         end

%         % Interpolates Distance to bedrock for the centre of each lake
%         BedrockDist = mapinterp(Bedrock_distance,C,[lake(:).centroidX],[lake(:).centroidY]);
%         for m = 1:length(BedrockDist)
%             lake(m).Rock_Distance = double(BedrockDist(m));
%         end
        
        % Assigns for each lake whether it is grounded (1) of on floating ice (0)
        grounded = inpolygon(cat(1,lake(:).centroidX),cat(1,lake(:).centroidY),Groundingline.X,Groundingline.Y);
        for n = 1:length(lake)
            lake(n).Grounded = double(grounded(n));
        end
        
    toc;
    
    shapewrite(lake,OutputFilePath);
    disp('-Saving shapefile...')
    % toc;
    
    else
        disp('No lakes are above the minimum area threshold, moving on to next time window...')
    end
    
    clearvars lake shapes new_output
    
    disp(['------------------------------------------------------------']);
    
end

else 
    disp('Not running Stage 4')
end

disp('Finished :-)');

warning('on')
%toc;
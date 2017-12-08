%This script counts cell nuclei per file
%BH 11/2015

close all
clear all
current_dir = cd;
%output file path
output_mean_path = strcat(current_dir,'\counter_output\cell_counter_mean.csv');
output_individual_path = strcat(current_dir,'\counter_output\cell_counter_');

%Cell-Count------------------------------------------------------------
prompt_cell_count = 'Do you wish to count cells?[y/n] ';
ans_cell_count = input(prompt_cell_count,'s');

if strcmp(ans_cell_count,'y') == 1

%prompts
prompt_show_fig = 'Do you wish to show figures?[y/n] ';
ans_show_fig = input(prompt_show_fig,'s');    
prompt_save = 'Do you wish to save output?[y/n] ';
ans_save = input(prompt_save,'s');

mkdir('counter_output');

%gets file in question
[FileName,PathName,FilterIndex] = uigetfile({'*.tif'},'MultiSelect','on');
cd(PathName); %changing directory to get files

%manual vs. automatic thresholding
prompt_region_search = 'Do you wish to manually select a region intensity for thresholding?[y/n] ';
ans_region_search = input(prompt_region_search,'s');

%number of files grabbed
switch class(FileName)
    case 'char'
        number_files = 1;
        file_flag = 1;
    case 'cell'
        number_files = numel(FileName);
        file_flag = 2;
end

%start nuclei indexing
previous_nuclei = 0;
%read images and save output
for i=1:number_files

    %indexing
    if file_flag == 1
        A = imread(FileName);
    elseif file_flag == 2
        A = imread(FileName{i});
    end
    B = imresize(A,0.5); %resize if too large

    if strcmp(ans_show_fig,'y') == 1
        figure, imshow(B);
        title('Raw file');
    end
    
    if strcmp(ans_region_search,'y') == 1
        %finding region where you want to zoom to grab colors
        rect = getrect(figure, imshow(B)); %allows you to select rectangular region [xmin ymin width height]
        disp(rect);
        rrect_1 = round(rect(1));
        rrect_2 = round(rect(2));
        rrect_3 = round(rect(3));
        rrect_4 = round(rect(4));

        B_rect = imcrop(B, [rrect_1 rrect_2 rrect_3 rrect_4]); %[xmin ymin width height]
        figure, imshow(B_rect);

        %grabbing ginput coordinates
        rgb_coord = ginput(1);
        x = rgb_coord(1);
        y = rgb_coord(2);
        xi = round(rect(1) + (x));
        yi = round(rect(2) + (y));
        disp(xi);
        disp(yi);

        %finding pixel intensities of ginput coordinates

        if size(B,3) == 3 %ginput intensity for a rgb image

            rgbValues = zeros(1,3);
            rgbValues(1,:) = reshape(B(yi(1), xi(1), :), 1, 3); %grabs array unit y,x,z from 1:3 and rearranges z elements in 1x3 array
            %code borrowed from http://iswwwup.com/t/0ac3d61e6e87/matlab-finding-saving-rgb-value-of-some-selected-pixels-of-an-image.html

            sprintf('The rgb values selected are: %d %d %d', rgbValues(1,1), rgbValues(1,2), rgbValues(1,3))

            %thresholding intensity for rgb values
            red_less = B(:,:,1) < rgbValues(1)+20; %value of red < red value of selected point
            red_more = B(:,:,1) > rgbValues(1)-20; %value of red > red value of selected point
            green_less = B(:,:,2) < rgbValues(2)+20;
            green_more = B(:,:,2) > rgbValues(2)-20;
            blue_less = B(:,:,3) < rgbValues(3)+20;
            blue_more = B(:,:,3) > rgbValues(3)-20;

            %2D array of w x h x 1 of colored pixels matching
            color_pixels = red_less & red_more & green_less & green_more & blue_less & blue_more;

            %converts a matrix to a grayscale image (0 = black, 1 = white)
            C = mat2gray(color_pixels);
            if strcmp(ans_show_fig,'y') == 1
                figure, imshow(C);
            end
            
        elseif size(B,3) == 1 %ginput intensity for a greyscale/BW image

            grayValues = zeros(1,1);
            grayValues(1,:) = reshape(B(yi(1), xi(1), 1), 1, 1); %grabs array unit y,x,z from 1:3 and rearranges z elements in 1x1 array

            sprintf('The grayscale value selected is: %d', grayValues)

            %thresholding intensity
            gray_less = B(:,:,1) < grayValues(1) + 80;
            gray_more = B(:,:,1) > grayValues(1) - 80;
            gray_pixels = gray_less & gray_more;

            C = im2bw(gray_pixels,0.3);
            if strcmp(ans_show_fig,'y') == 1
                figure, imshow(C);
            end
        end

    elseif strcmp(ans_region_search,'n') == 1 %automatically select an intensity for thresholding based on mean value

        if size(B,3) == 3 %ginput intensity for a rgb image
          
            BB = im2double(B);
            BB_red = reshape(BB(:,:,1),1,numel(BB(:,:,1))); %reshapes matrix to 1xm vector
            BB_green = reshape(BB(:,:,2),1,numel(BB(:,:,2))); %reshapes matrix to 1xm vector
            BB_blue = reshape(BB(:,:,3),1,numel(BB(:,:,3))); %reshapes matrix to 1xm vector
            mean_intensity_red = mean(BB_red)*2550;
            mean_intensity_green = mean(BB_green)*2550;
            mean_intensity_blue = mean(BB_blue)*2550;

            sprintf('The rgb values selected are: %d %d %d', mean_intensity_red, mean_intensity_green, mean_intensity_blue)

            %thresholding intensity
            red_less = B(:,:,1) < mean_intensity_red+20; %value of red < red value of selected point
            red_more = B(:,:,1) > mean_intensity_red-20; %value of red > red value of selected point
            green_less = B(:,:,2) < mean_intensity_green+20;
            green_more = B(:,:,2) > mean_intensity_green-20;
            blue_less = B(:,:,3) < mean_intensity_blue+20;
            blue_more = B(:,:,3) > mean_intensity_blue-20;
            %2D array of w x h x 1 of colored pixels matching
            color_pixels = red_less & red_more & green_less & green_more & blue_less & blue_more;

            %converts a matrix to a grayscale image (0 = black, 1 = white)
            C = mat2gray(color_pixels);
            if strcmp(ans_show_fig,'y') == 1
                figure, imshow(C);
            end
            
        elseif size(B,3) == 1 %ginput intensity for a greyscale/BW image
            
            %BB_adj = adapthisteq(B);
            %BB = im2double(BB_adj);
            BB = im2double(B);
            BB_gray = reshape(BB,1,numel(BB));
            mean_intensity_gray = mean(BB_gray)*2550;
            min_intensity_gray = min(BB_gray)*2550;
            max_intensity_gray = max(BB_gray)*2550;
            
            sprintf('The grayscale value selected is: %d', mean_intensity_gray)

            %thresholding intensity
            %gray_less = B(:,:,1) < mean_intensity_gray + abs(mean_intensity_gray-max_intensity_gray)/2;
            %gray_more = B(:,:,1) > mean_intensity_gray - abs(mean_intensity_gray-min_intensity_gray)/8;
            %gray_pixels = gray_less & gray_more;
            
            %gray_raw = B(:,:,1) > mean_intensity_gray;
            
            %Z = zscore(BB);
            Z = zscore(BB);
                        
            gray_pixels = Z > 3;
            
            %converts a matrix to a grayscale image (0 = black, 1 = white)
            C = im2bw(gray_pixels, 0);
            if strcmp(ans_show_fig,'y') == 1
                figure, imshow(C);
            end
            
        end

    end

    %dilation and filling holes
    %D = bwmorph(C, 'dilate');
    
    %cd(current_dir)
    %C_bpass=  bpass(B,2,10);
    %cd(PathName);
    %D = edge(C_bpass,'canny');
    C = bwareaopen(C, 50);
    %D = edge(C,'canny',[mean_intensity_gray/(max_intensity_gray+1) (max_intensity_gray/(max_intensity_gray+1))]);
    D = edge(C,'canny');
    %D = edge(C,'sobel');
    %D = C;
    E = imfill(D, 'holes');
    F = bwareaopen(E, 50); %threshold out anything <50 pixels
    
    %E = imfill(bwmorph(bwmorph(bwmorph(D, 'spur'), 'dilate'), 'spur'), 'holes');
    
    %gives stats about cells identified
    stats = regionprops(F, 'All');

    %premake matrices
    nuclei_area = zeros(numel(stats),1);
    nuclei_orientation = zeros(numel(stats),1);
    nuclei_eccentricity = zeros(numel(stats),1);
    centroid_x = zeros(numel(stats),1);
    centroid_y = zeros(numel(stats),1);
        
    if strcmp(ans_show_fig,'y') == 1
        subplot(2,2,1), imshow(D);
        title('edge');
        subplot(2,2,2), imshow(E);
        title('edge, fill holes');
        subplot(2,2,3), imshow(F);
        title('edge, fill holes, remove objects');
    end
    
    %inputs stats data into premade matrices
    for j=1:numel(stats)
        %disp(stats(j));
        nuclei_area(j)=stats(j).Area;
        nuclei_orientation(j) = stats(j).Orientation;
        nuclei_eccentricity(j) = stats(j).Eccentricity;
        centroid_x(j) = stats(j).Centroid(1);
        centroid_y(j) = stats(j).Centroid(2);
        
        %nuclei index (nth nuclei found)
        nuclei_index = j + previous_nuclei;
        
        %saving individual image outputs to file...
        image_num = num2str(i); %processed imaged number
        current_img_path = strcat(output_individual_path,image_num,'.csv');
        if strcmp(ans_save,'y') == 1
            if exist(current_img_path, 'file') == 0 %does not exist
                fid_indiv = fopen(current_img_path, 'w+');
            elseif exist(current_img_path, 'file') == 2 %exists
                fid_indiv = fopen(current_img_path, 'a+');
            end
            fprintf(fid_indiv,'%s%s%s%s%s%s', num2str(nuclei_index),',', num2str(centroid_x(j)),',', num2str(centroid_y(j)),',', num2str(nuclei_area(j)),',', num2str(nuclei_orientation(j)),',', num2str(nuclei_eccentricity(j)));
            fprintf(fid_indiv,'\n');
            fclose(fid_indiv);
        end
    end
    clearvars j current_img_path;
    
    if strcmp(ans_show_fig,'y') == 1
        %plot centroids
        figure, imshow(F);
        hold on
        for j=1:numel(stats)
            %plot(imgca, centroid_x(j), centroid_y(j), 'r*')
            %label number
            
            text(centroid_x(j), centroid_y(j), sprintf('%d', j+previous_nuclei), 'Color', 'r','FontWeight', 'bold', 'FontSize', 8);
        end
    end
    
    sprintf('The number of nuclei found is: %d', numel(stats))
    sprintf('The mean area of the nuclei is: %d', mean(nuclei_area))
    sprintf('The mean orientation of the nuclei is: %d', mean(nuclei_orientation))
    sprintf('The mean eccentricity of the nuclei is: %d', mean(nuclei_eccentricity))
    
    %saving mean output per image to file...
    if strcmp(ans_save,'y') == 1
        if exist(output_mean_path, 'file') == 0 %does not exist
            fid_mean = fopen(output_mean_path, 'w+');
            fprintf(fid_mean,'%s%s%s', num2str(numel(stats)),',', num2str(mean(nuclei_orientation)),',', num2str(mean(nuclei_eccentricity)));
            fprintf(fid_mean,'\n');
            fclose(fid_mean);
        elseif exist(output_mean_path, 'file') == 2 %exists
            fid_mean = fopen(output_mean_path, 'a+');
            fprintf(fid_mean,'%s%s%s', num2str(numel(stats)),',', num2str(mean(nuclei_orientation)),',', num2str(mean(nuclei_eccentricity)));
            fprintf(fid_mean,'\n');
            fclose(fid_mean);
        end
    end
    
    %save image
    imwrite(F,strcat(output_individual_path,'img_',num2str(i),'.tif'),'tif');
    
    %nuclei index
    previous_nuclei = previous_nuclei + numel(stats);
    
end
clearvars i number_files file_flag centroid_x centroid_y;

%total nuclei found
sprintf('The total number of nuclei (inc. degenerate) found is: %d', previous_nuclei)
clearvars previous_nuclei;

end
%Cell-Count-End---------------------------------------------------------

%Cell-Depth-Analyis-----------------------------------------------------
cd(current_dir);
prompt_analysis = 'Do you wish to average each cell depth?[y/n] ';
ans_analysis = input(prompt_analysis,'s');

if strcmp(ans_analysis,'y') == 1
        
    %gets file in question
    [FileName,PathName,FilterIndex] = uigetfile({'*.csv'},'MultiSelect','on');
    cd(PathName); %changing directory to get files
    %depth and step size prompts
    prompt_step_size = 'What is your step size (x-axis [um])? ';
    ans_step_size = input(prompt_step_size);
    prompt_start_depth = 'What is your starting depth [um]? ';
    ans_start_depth = input(prompt_start_depth);
    
    %number of centroid files selected
    switch class(FileName)
        case 'char'
            number_files = 1;
            file_flag = 1;
        case 'cell'
            number_files = numel(FileName);
            file_flag = 2;
    end
    
    %grabbing centroid files
    csv_centroid = cell(number_files,1);
    for i=1:number_files
        if file_flag == 1 %single centroid files grabbed
            fid_raw_cell = fopen(FileName,'r');
            csv_centroid{i,1} = textscan(fid_raw_cell,'%d%d%d%d%d%d','Delimiter',',');
            fclose(fid_raw_cell);
        elseif file_flag == 2 %multiple centroid files
            fid_raw_cell = fopen(FileName{i},'r');
            csv_centroid{i,1} = textscan(fid_raw_cell,'%d%d%d%d%d%d','Delimiter',',');
            fclose(fid_raw_cell);
        end
    end
    clearvars i;
        
    %premake empty tracking array (double)
    tracking_array_raw = [];
    tracking_array = cell(0);
    previous_nuclei = 0;
    
    %grab directory, count files --> for numel(files) take i, i+1 for
    %comparison
    
    %matching up centroids of previous n=2 images for all frames
    for i=1:number_files
        %x,y centroid to be compared
        for j=1:numel(csv_centroid{i}{1})
            
            centroid_x = double(csv_centroid{i}{2}(j));
            centroid_y = double(csv_centroid{i}{3}(j));
            next_centroid = i+1;
            nuclei_index = j + previous_nuclei;
            
            if next_centroid < number_files + 1 %when to stop
                
                nuclei_matched = 0; %counter for nuclei matching
                for k=1:numel(csv_centroid{next_centroid}{1})
                    
                    disp('Calculating distance...');
                    distance_x = double(csv_centroid{next_centroid}{2}(k))-centroid_x;
                    distance_y = double(csv_centroid{next_centroid}{3}(k))-centroid_y;
                    distance = sqrt(distance_x^2 + distance_y^2);
                    disp(distance);
                    
                    if distance < 2 %<2 pixel difference = match
                        sprintf('nuclei %d matches with %d', csv_centroid{i}{1}(j), csv_centroid{next_centroid}{1}(k))
                        tracking_array_raw(nuclei_index, 1) = csv_centroid{i}{1}(j);
                        tracking_array_raw(nuclei_index, 2) = csv_centroid{next_centroid}{1}(k);
                        nuclei_matched = nuclei_matched + 1;
                    end
                    %calculate distance here
                    %be wary of non-existent indices
                end
                if nuclei_matched == 0 %if no nuclei were matched
                    tracking_array_raw(nuclei_index, 1) = csv_centroid{i}{1}(j);
                    tracking_array_raw(nuclei_index, 2) = 0;
                end
                    
            end
            sprintf('x: %d, y: %d', centroid_x, centroid_y)
        end
        
        previous_nuclei = previous_nuclei + numel(csv_centroid{i}{1});
        
    end
    clearvars i j k previous_nuclei;
    
    %organizing tracking array
    for i=1:numel(tracking_array_raw(:,1))
        
        track_key = tracking_array_raw(i,2);
        %creates tracking array, seeding first elements with initial vector
        tracking_array{i}(1:2) = tracking_array_raw(i,:);
        previous_tracking = numel(tracking_array_raw(i,:)) + 1;
        %previous_tracking = 0;
        
        if track_key ~=0 %remove searching for non-existant indices (0)
            for j=i+1:numel(tracking_array_raw(:,1))
                
                 if ismember(track_key, tracking_array_raw(j,:)) == 1
                     index_tracking = previous_tracking + 1;
                     
                     %appending to array
                     tracking_array{i}(index_tracking:index_tracking-1+numel(tracking_array_raw(j,:))) = tracking_array_raw(j,:);
                     
                     if tracking_array_raw(j,2) ~= 0 %exclude searching for non-existant indices (0)
                        track_key = tracking_array_raw(j,2); %track key is constantly changing
                     end
                     
                     %clearing (zeroing) raw array so that it doesn't need
                     %to be used later
                     tracking_array_raw(j,:) = 0;
                     
                     %update tracking array size
                     previous_tracking = numel(tracking_array{i});    
                 end
                 
            end
        end
        
    end
    clearvars i j previous_tracking;
    
    %depth tracking
    start_depth = ans_start_depth; %depth, specified
    step_size = ans_step_size; %step size, specified
    depth_array = cell(number_files,1);
    %creating depth reference array to calc. mean depth of nuclei
    for i=1:number_files
        depth_array{i,1}(1) = min(csv_centroid{i}{1}); %depth_array cell min
        depth_array{i,1}(2) = max(csv_centroid{i}{1}); %depth_array cell max
    end
    clearvars i prompt_step_size prompt_start_depth ans_step_size ans_start_depth;
    
    %premake arrays
    tracking_array_depth = cell(0);
    tracking_array_depth_mean = cell(0);
    tracking_array_str = cell(0);
    tracking_array_depth_str = cell(0);
    
    %nuclei index depth conversion
    for i=1:numel(tracking_array)
        %cleaning up tracking array
        tracking_array{i} = unique(tracking_array{i}); %unique values
        tracking_array{i}(tracking_array{i}==0) = []; %removing zeros
        %converting index values to depth values        
        for j=1:numel(tracking_array{i})
            for k=1:number_files
                %if tracking array nuclei number between range
                if tracking_array{i}(j)>=depth_array{k,1}(1) && tracking_array{i}(j)<=depth_array{k,1}(2)
                    tracking_array_depth{i}(j) = ((k-1)*step_size) + start_depth; %convert to depth
                end
            end
        end
    end
    clearvars i;
    
    %calculating individual depth by averaging depth values and saving output
    current_depth_path_raw = strcat(output_individual_path,'depth_raw.csv');
    for i=1:numel(tracking_array_depth)
        %mean of depth
        tracking_array_depth_mean{i} = round(mean(tracking_array_depth{i})); %round to nearest um
        sprintf('Mean depth: %d', tracking_array_depth_mean{i})
        
        tracking_array_depth_mean_dbl(i) = tracking_array_depth_mean{i};
        
        %string for saving raw frames alongside depth
        tracking_array_str{i} = strrep(num2str(tracking_array{i}),' ','-');
        tracking_array_depth_str{i} = strrep(num2str(tracking_array_depth{i}),' ','-');
        
        %saving
        if exist(current_depth_path_raw, 'file') == 0 %does not exist
            fid_depth = fopen(current_depth_path_raw, 'w+');
            fprintf(fid_depth,'%s%s%s', num2str(tracking_array_depth_mean{i}),',', tracking_array_str{i},',', tracking_array_depth_str{i});
            fprintf(fid_depth,'\n');
            fclose(fid_depth);
        elseif exist(current_depth_path_raw, 'file') == 2 %exists
            fid_depth = fopen(current_depth_path_raw, 'a+');
            fprintf(fid_depth,'%s%s%s', num2str(tracking_array_depth_mean{i}),',', tracking_array_str{i},',', tracking_array_depth_str{i});
            fprintf(fid_depth,'\n');
            fclose(fid_depth);
        end
        
    end
    clearvars i;
    
    %summing total depth and saving output
    current_depth_path_mean = strcat(output_individual_path,'depth_mean.csv'); 
    tracking_array_depth_mean_dblu = unique(tracking_array_depth_mean_dbl); %unique values
    tracking_array_depth_mean_dblu(isnan(tracking_array_depth_mean_dblu)) = []; %removing NaN values
    
    for i=start_depth:start_depth+number_files-1 %specifying starting/ending depth
        if ismember(i*step_size, tracking_array_depth_mean_dblu) == 1
            depth_number = sum(ismember(tracking_array_depth_mean_dbl, i*step_size)); %summing occurences
        elseif ismember(i*step_size, tracking_array_depth_mean_dblu) == 0
            depth_number = 0; %if no nuclei are tracked to given depth
        end
        
        sprintf('Depth: %d Nuclei: %d', i*step_size, depth_number)
        
        if exist(current_depth_path_mean, 'file') == 0 %does not exist
            fid_depth_mean = fopen(current_depth_path_mean, 'w+');
            fprintf(fid_depth_mean,'%s', num2str(depth_number));
            fprintf(fid_depth_mean,'\n');
            fclose(fid_depth_mean);
        elseif exist(current_depth_path_mean, 'file') == 2 %exists
            fid_depth_mean = fopen(current_depth_path_mean, 'a+');
            fprintf(fid_depth_mean,'%s', num2str(depth_number));
            fprintf(fid_depth_mean,'\n');
            fclose(fid_depth_mean);
        end     

    end
    clearvars i;

end
%Cell-Depth-Analyis-End-------------------------------------------------

%Pericyte-Depth---------------------------------------------------------
cd(current_dir);
prompt_pericyte = 'Do you wish to threshold pericytes?[y/n] ';
ans_pericyte = input(prompt_pericyte,'s');

if strcmp(ans_pericyte,'y') == 1
    %calculate starting depth of pericytes via VE-Cadherin staining
    cd('C:\Users\Admin\Desktop\Matlab_scripts\nico\examples\flow\file_1\analysis');
    %gets file in question
    [FileName,PathName,FilterIndex] = uigetfile({'*.tif'},'MultiSelect','on');
    cd(PathName); %changing directory to get files
    %number of files grabbed
    switch class(FileName)
        case 'char'
            number_files = 1;
            file_flag = 1;
        case 'cell'
            number_files = numel(FileName);
            file_flag = 2;
    end

    %start nuclei indexing
    previous_nuclei = 0;
    %read images and save output
    for i=1:number_files

        %indexing
        if file_flag == 1
            A = imread(FileName);
        elseif file_flag == 2
            A = imread(FileName{i});
        end
        B = imresize(A,0.5); %resize if too large
        C = mat2gray(B);
        D = reshape(C,1,numel(C));
        sprintf('The sum of pixel intensities = %d', sum(D))
        intensities(i) = sum(D);
        
    end
    
    %crop first few microns
    [min_int index_min] = min(intensities(10:60));
    index_min = index_min+10;
    
    x_intensities = 1:numel(intensities);
    figure, plot(x_intensities,intensities(x_intensities));
    xlabel('depth [um]','FontSize',16);
    ylabel('intensity [-]','FontSize',16);
    hold on
    text(0, 500, sprintf('Threshold at: %d um', index_min), 'Color', 'r','FontWeight', 'bold', 'FontSize', 12);

end
%Pericyte-Depth-End-----------------------------------------------------


%Cell-Depth-Graph-------------------------------------------------------
cd(current_dir);
prompt_graph = 'Do you wish to graph depth vs cell number?[y/n] ';
ans_graph = input(prompt_graph,'s');

if strcmp(ans_graph,'y') == 1
    disp('Select your processed file (raw or averaged)');
    
    [FileName_counter,PathName_counter,FilterIndex] = uigetfile('*.csv');
    path = strcat(PathName_counter,'\',FileName_counter);
    cd(PathName_counter);
    
    fid_counter = fopen(path, 'r');
    %csv_counter = textscan(fid_counter,'%d%d%d','Delimiter',','); %for each column, need formatSpec
    csv_counter = textscan(fid_counter,'%d','Delimiter',','); %for each column, need formatSpec
    fclose(fid_counter);
    
    prompt_step_size = 'What is your step size (x-axis) [um]? ';
    ans_step_size = input(prompt_step_size);
    
    %looking at truncated section (e.g. 1 to 100 um)
    disp(['Min Section: ', num2str(ans_step_size), ' Max Section: ', num2str(numel(csv_counter{1})*ans_step_size)]);
    prompt_section = 'What section [um] would you like to examine (separate min/max with comma or enter max)? ';
    ans_section = input(prompt_section,'s');
    
    ans_section = strrep(ans_section,',',' ');
    ans_section = str2num(ans_section);
    
    if numel(ans_section) == 2
        if strcmp(ans_pericyte,'y') == 1
            section_min = index_min;
        else
            section_min = ans_section(1);
        end
        section_max = ans_section(2);
    elseif numel(ans_section) == 1
        if strcmp(ans_pericyte,'y') == 1
            section_min = index_min;
        else
            section_min = 1; %default 1 um starting
        end
        section_max = ans_section(1);
    elseif numel(ans_section) > 2
        error('Unable to identify section');
    end

    %convert raw data to proportion
    array = double(csv_counter{1});
    array_sum = sum(array);
    array_norm = array/array_sum;
    
    %setting x-axis length, values must be integers
    depth = linspace(1, numel(csv_counter{1})*ans_step_size, numel(csv_counter{1}));
    depth_truncated = depth((section_min/ans_step_size):(section_max/ans_step_size));
    array_truncated = array_norm((section_min/ans_step_size):(section_max/ans_step_size));

    [max_array index] = max(array_norm);
    sprintf('The maximum number of cells is at %d um', index)

    %graph
    bar(depth_truncated, array_truncated);
    %axis([minx maxx miny maxy]);
    axis([section_min max(depth_truncated) min(array_norm) max(array_norm)]);
    set(gca,'FontSize', 24);
    xlabel('depth [um]','FontSize',16);
    ylabel('% of cells (nuclei)','FontSize',16);
end
%Cell-Depth-Graph-End---------------------------------------------------

cd(current_dir);
disp('End script...goodbye');
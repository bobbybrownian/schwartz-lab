%This script splits tiff images and outputs files in a new folder
%BH 9/2015, v1.1 production

%--INSTRUCTIONS---------------------------------------------------------

%--INSTRUCTIONS-END-----------------------------------------------------

close all
clear all

current_dir = cd;

%--SPLIT-TIFF-PROMPT----------------------------------------------------
%prompts to split tiff and create directories, if 'n', then proceed to
%processing movieData.mat
prompt_split_tiff = 'Would you like to split a tiff? [y/n] ';
split_tiff = input(prompt_split_tiff,'s');
switch split_tiff
    
%if split_tiff == 'y'
case 'y'   
%--SPLIT-TIFF-PROMPT-END------------------------------------------------

%--SETTING-DIRECTORY----------------------------------------------------
%creating directory, would need to do once per session
disp('Select a directory to work in');
dir_path = uigetdir;

prompt_name_dir = 'What would you like to name your directory? ';
directory_name = input(prompt_name_dir,'s');

mkdir(dir_path,directory_name); %path, name
%--SETTING-DIRECTORY-END------------------------------------------------


%--SETTING-CHANNELS-----------------------------------------------------
prompt_num_channel = 'How many channels would you like to analyze? ';
channel_num = input(prompt_num_channel);

channel_array = {channel_num,1,1}; %3 rows, 1 column, 1 z

for n=1:channel_num
    n_str = num2str(n);
    prompt_name_channel = ['Name of channel ',n_str,'? ']; %strcat removes white spaces
    channel_array{n} = input(prompt_name_channel,'s');
end
clearvars n

%creating cell subfolders
prompt_cell_number = 'How many files are you analyzing? ';
prompt_folder_naming = 'How would you like to name your folders? ';
cell_folder_num = input(prompt_cell_number);
folder_naming = input(prompt_folder_naming,'s');

for i=1:cell_folder_num
    i_var = num2str(i);
    cell_folder_name = strcat(folder_naming,'_',i_var);
    sub1folder_path = strcat(dir_path,'\',directory_name);
    mkdir(sub1folder_path,cell_folder_name);
    sub2folder_path = strcat(dir_path,'\',directory_name,'\',cell_folder_name);
    for n=1:channel_num
        mkdir(sub2folder_path,channel_array{n});
        disp(channel_array{n});
    end    
    disp(cell_folder_name);
end
clearvars n;

disp('Directories created...');
dir_path_cell = strcat(dir_path,'\',directory_name);
cell_folder_class = folder_naming;
%--SETTING-CHANNELS-END-------------------------------------------------


%--SELECTING-TIFF-FILES-------------------------------------------------
%need to split tiff stack
disp('Select tiff file(s) to split');

%[FileName,PathName,FilterIndex] = uigetfile({'*.tif'});
[FileName,PathName,FilterIndex] = uigetfile({'*.tif'},'MultiSelect','on');

cd(PathName); %cd to file location

%identifying number of tifs selected...
switch class(FileName)
    case 'char'  %if FileName is character array --> only 1 tif selected
        num_files = 1;
        info = imfinfo(FileName);
        num_images{1} = numel(info);
        tif_multiple = 0;
        %saving information text file
        info_path = strcat(dir_path_cell,'\',cell_folder_class,'_1','\info.txt');
        fid_info = fopen(info_path,'w+'); %creates new writable file      
        fprintf(fid_info,'%s\r\n','Original file name:',FileName,' ','Channels:');%fileID,format,data1,data2,... (w/newline)
            for n=1:channel_num
                fprintf(fid_info,'%s\r\n',channel_array{n});%fileID,format,data1,data2,... (w/newline)
            end
        fclose(fid_info);
    case 'cell' %if FileName is cell array --> 2+ tif selected
        num_files = numel(FileName);
        info = {num_files};
        num_images = {num_files};
        tif_multiple = 1;
        %tabulating files selected
        for f=1:num_files
            f_var = num2str(f);
            info{f} = imfinfo(FileName{f});
            num_images{f} = numel(info{f});
            %disp(info{f});
            %saving information text file
            info_path = strcat(dir_path_cell,'\',cell_folder_class,'_',f_var,'\info.txt');
            fid_info = fopen(info_path,'w+'); %creates new writable file
            fprintf(fid_info,'%s\r\n','Original file name: ',FileName{f},' ','Channels:');%fileID,format,data1,data2,... (w/newline)
            for n=1:channel_num
                fprintf(fid_info,'%s\r\n',channel_array{n});%fileID,format,data1,data2,... (w/newline)
            end
            fclose(fid_info);
        end
end

%if number of files selected not equal to number of cells inputted
if cell_folder_num ~= num_files
    error('The number of files selected != number cells analyzing');
end
%--SELECTING-TIFF-FILES-END---------------------------------------------


%--SPLITTING-AND-SAVING-TIFF-STACKS-------------------------------------
for t=1:num_files
%----------------

%getting folder path for cell being analyzed
t_var = num2str(t);
cell_folder_path  = strcat(dir_path_cell,'\',cell_folder_class,'_',t_var);

%splitting and saving all the images
for j=1:(num_images{t}/channel_num)
    
    z = 1; %counter to reiterate every frame (block)
    
    for k=(((j-1)*channel_num)+1):(((j-1)*channel_num)+(channel_num))
        %seeing if image stack is part of a set (if more than 1 tif
        %selected)
        switch tif_multiple
            case 1
            A = imread(FileName{t}, k, 'Info', info{t});
            case 0
            A = imread(FileName, k, 'Info', info);
        end
        
        k_var = num2str(k);
        sprintf('Image %d of %d', k, num_images{t})

        imshow(A);
        %setting 1000x1000 image to save (without borders)
        set(gca,'position',[0 0 1 1],'units','normalized');
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10]);            
        index_channel = ((k-z)/channel_num)+1;
        %adding zeros to name of string
        if index_channel<10
            prepend_zeros = '000';
        else
            if index_channel>=10
                if index_channel<100
                    prepend_zeros = '00';
                end
                if index_channel>=100
                    prepend_zeros = '0';
                end
            end
        end
        index_channel = num2str(index_channel);
        cell_folder_path_channel = strcat(cell_folder_path,'\',channel_array{z},'\',channel_array{z},prepend_zeros,index_channel,'.tif');
        imwrite(A,cell_folder_path_channel,'tif');
        disp([channel_array{z} ' splitting...']);            
        
    %disp(z);
    z = z+1;
    end
    
    disp('end frame');
    
end

clearvars cell_folder_path;
%----------------
end
%--SPLITTING-AND-SAVING-TIFF-STACKS-END---------------------------------
 
end
%--SPLIT-TIFF-PROMPT-END------------------------------------------------

disp('End script...goodbye');
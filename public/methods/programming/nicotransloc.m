function tres = nicotransloc(stn, nsn)
% Nicotransloc quantifies the translocation of a stain into the nucleus.
% The TF (translocation factor) is a unitless value that is normalized such that
% a cell with total translocation has a value of one while a cell without any translocation
% has a value of zero.
%
% The translocation factor is calculated by taking the product of the both the area
% and the intensity of the stain present in the nucleus and dividing it by
% the product of the total stain area and the total intensity of the stain.
%
% INPUT
% stn - file name of tiff stack containing stain of interest
% nsn - file name of tiff stack containing nuclear stain
% OUTPUT
% tres - A cell containing three different entries:
%
%        1. A list of all the raw data for each individual cell with
%        columns as follows: [Stain Area, Total Stain Intensity in
%        Nucleus, Stain Area within Nucleus, Total Stain Intensity].
%
%        2. The mean translocation factor of the stack.
%
%        3. The standard error of the translocation factor of the stack.
%
%
% Example Input: syf = nicotransloc('p65*syn*ow.tif','dapi*syn*ow.tif');
% 
% 
% 7-27-2012
% Writen by Tyler Ross 
% tyler.ross@yale.edu
% Martin A. Schwartz Lab 


% Find files
sfn = dir(stn);
nfn = dir(nsn);

% Find number of pages in stack
ns = numel(imfinfo(sfn(1).name));
nn = numel(imfinfo(nfn(1).name));
ares=[];
if ns ~= nn,
    error('Number of images do not match!');
end

% Loop through pages of stack
for i = 1:ns,
    
    % Read in images
    sim = double(imread(sfn(1).name,i));
    nim = double(imread(nfn(1).name,i));
    
    % Normalize Images
    nsim = sim/max(sim(:));
    nnim = nim/max(nim(:));

    %Apply filter to images to differentiate signal from background
    nsim = nsim.*imclose(imopen((nsim>mean(nsim(:))),ones(5,5)),ones(5,5));
    nnim = nnim.*imfill(imopen(imclose(nnim>graythresh(nnim),ones(5,5)),ones(5,5)),'holes');
    
    
    %Remove small objects in both images
    smsk = bwlabel(nsim,4);
    props = regionprops(smsk,'Area');
    szs = cat(1,props.Area);
    bsz = find(szs > 750);
    nnmsk = smsk*0;
    for j = 1:numel(bsz),
        nnmsk(smsk==bsz(j))=1;
    end
    smsk = nnmsk;
    smsk = bwlabel(smsk,4);
    nsim = nsim.*(smsk>0);
    
    nmsk = bwlabel(nnim,4);
    props = regionprops(nmsk,'Area');
    szs = cat(1, props.Area);
    bsz = find(szs > 500);
    nnmsk = nmsk*0;
    for j = 1:numel(bsz),
        nnmsk(nmsk==bsz(j))=1;
    end
    nmsk = nnmsk;
    nmsk = bwlabel(nmsk,4);
    
    
    %Eliminate stains that are not associated with a nucleus
    cmsk = bwlabel(imfill(smsk,'holes'),4);
    for j=1:max(cmsk(:)),
        n = (cmsk==j).*nmsk;
        if max(n(:)) == 0,
            cmsk(cmsk==j) = 0;
        end
        
    end
    nsim = nsim.*(cmsk>0);
    
    
    %Segment cells and remove ones that are on the image boarder
    smsk = bwlabel(imclearborder(fseg(nsim,750,nmsk)).*(smsk>0),4);
    
    %Link IDs between nucleus and translocating stain
    nmsk = nmsk.*-1;
    for j=1:max(abs(nmsk(:))),
        n = (nmsk == -1*j).*smsk;
        if max(n(:)) >0,
            nmsk(nmsk(:)== -1*j) = max(n(:));
        else
            nmsk(nmsk(:) == -1*j) = 0;
        end
        
    end
    
    
    %Mask of stain inside nucleus
    imsk = (nmsk >0).*smsk;

    %Quantify Area and Total Intensity for TF
    props = regionprops(smsk,'Area');
    propst = regionprops(imsk,'Area');
    u = unique(nmsk);
    u=u(2:end);
    for j =1:numel(u),
        res{j,1} = props(u(j)).Area;
        res{j,2} = sum(sim(imsk == u(j)));
        res{j,3} = propst(u(j)).Area;
        res{j,4} = sum(sim(smsk == u(j)));
    end
    res = cell2mat(res);
    ares = [ares;res];
    clear('res');
end
% Calculate TF
iv = (ares(:,2).*ares(:,3));
ev = (ares(:,4).*ares(:,1));
tres{1,1} = ares;
tres{1,2} = mean(iv./(ev));
tres{1,3} = std(iv./(ev))/sqrt(size(ares,1));

end
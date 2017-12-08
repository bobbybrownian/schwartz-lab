close all;
clear all;
img_direc_base = 'D:\Keiichiro\2017-01-25 FRET Gi2 serum paired & timecourse\Analysis (488nm_25per,561nm_2per)\';
gfpgrad=imread(strcat(img_direc_base,'Atto488 60x- grad.tif'));
rfpgrad=imread(strcat(img_direc_base,'Atto565 60x- grad.tif'));

%fretgrad=imread(strcat(img_direc_base,'Atto488 100x1x_FRETgrad.tif'));
fretgrad=gfpgrad;
cy5flag=0; %=1 to get intesity of cy5 at FRET pixel else 0
if cy5flag==1
    cy5grad=imread(strcat(img_direc_base,'Atto565 100x-grad.tif'));
end
flag=0;
%0- if thresholded image present
%1-for manually cropping fret image for cytoplasmic pool or ROI
%2 if no cropping, full image required
%5- for if thresholded image present and need full nonFa region
if flag==5
    nonFafilepre='-nonFA';
else
    nonFafilepre='';
end
imgsizex=1000;
imgsizey=1000;
roino=10;   %make 10 for cytosol

subdirec='Gi2-FRET 60xlive timecourse3 before-';
subdirec1=subdirec;
%subdirec1=[subdirec,' ',sprintf('%0d',12)];

opfilepre='';
subdirecno=1;
%timepoints=6;
%Spinning disk 100x 
        %donleak- 0.05463; acccross-0.07617; 
%after 2016-02-15 spinning disk- 
        % donleak-0.0518; acccross- 0.03563
        % donleak=0.05312; acccross=0.03992;
%Spinning disk 60x 12% 488, 8% 565
        % donleak=0.05445; acccross=0.04041;
%Spinning disk 2017-02-15 100x 12% 488, 8% 565
        % donleak=0.05361; acccross=0.05158;
%Spinning disk 2017-02-16 100x 12% 488, 8% 565   after clean up
        % donleak=0.0517; acccross=0.05314;
%Spinning disk 100x 2X2 binning 12% 488, 8% 565-  
        % donleak=0.06916; acccross=0.06125;
%Ring TIRF 100x1x 
        % donleak- 0.08627; acccross-0.21189;
%Ring TIRF 60x1x?? 
        % donleak- 0.07758; acccross-0.11867;
%POLscope 
        % donleak=0.18009;    acccross=0.03381;
        % donleak=0.05445; acccross=0.04041;
%Spinning dish 12-2016 15% 488 1%RFP 
        % donleak=0.05836; acccross=0.3856;
 %Spinning dish 12-2016 25% 488 2%RFP 
        donleak=0.04956; acccross=0.33225;

smoothtype = fspecial('average'); %for 3x3 image smoothening
cmap=colormap('jet');
cmap(1,3)=0.0;
fretscalemin=0.00;
fretscalemax=0.2;

binmin=0.0; binmax=0.5; binsize=0.005;
edges=[binmin:binsize:binmax];
edgereal=[binsize/2:binsize:(binmax-binsize/2)];
edgereal=edgereal';

% Get the background intensity 
%img_direc = strcat(img_direc_base,subdirec1,'\');
%cd(img_direc);
%close all;
%time=1;
%time=strcat('T',sprintf('%05d',time));
%time=strcat('img_',sprintf('%09d',time));
               
%donfname = strcat(time,'_488-Laser-Shutter-Andor sCMOS Camera_000.tif');
%donfname = strcat(time,'_GFP (FRET)_000.tif');
%donfname = strcat(time,'C01Z001.tif');

%don=imread(donfname);
%don=imfilter(don,smoothtype);

%bkimg=imcrop(don, [0 500]);   %120 forRingTIRF
%avgbk=mean(mean(bkimg))

avgbk=450;%450;%465;460 

th=1.0;
intth=avgbk/10; %for ring tirf-avgbk/15; for EM
checkintth2=1;  %for dual camera no threshld image but just intensity threshold ...1 for extra int threshold else 0
intth2=avgbk/3;   %use avgbk/2; for stringent thresholding use avgbk; for ring tirf use avgbk/4;
intthmax= 12000.0;
objective=60;  % choose objective 60/100
if objective==60
    rowshift=0;
    colshift=0;
end
if objective==100
    rowshift=0;  %0     % CHANGE for DUAAL CAMERA RING TIRF and other microscope
    colshift=0;  %1
end
for sub=1:subdirecno    
    img_direc = strcat(img_direc_base,subdirec1,'\');
    cd(img_direc);    
    %timepoints=size(dir('img*_488-Laser-Shutter-Andor sCMOS Camera_000.tif'),1);
    timepoints=size(dir('T*C01Z001.tif'),1);
    %timepoints=1;    
    for i=1:timepoints        
        close all;
        time=i;  %ring tirf-i-1;
        time=strcat('T',sprintf('%05d',time));
        %time=strcat('img_',sprintf('%09d',time));               
        
        %donfname = strcat(time,'_488-Laser-Shutter-Andor sCMOS Camera_000.tif');
        %accfname = strcat(time,'_488-Laser-Shutter-Andor sCMOS Camera-1_000.tif');
        %rfpfname = strcat(time,'_568-Laser-Shutter-Andor sCMOS Camera-1_000.tif');        
        
        %donfname = strcat(time,'_GFP (FRET)_000.tif');
        %accfname = strcat(time,'_GFP-RFP FRET_000.tif');
        %rfpfname = strcat(time,'_TagRFP_000.tif');        
        
        donfname = strcat(time,'C01Z001.tif');
        accfname = strcat(time,'C02Z001.tif');
        rfpfname = strcat(time,'C03Z001.tif');    
        
        don=imread(donfname);
        acc=imread(accfname);
        rfp=imread(rfpfname);
        
        don=double(don);
        acc=double(acc);
        rfp=double(rfp);
        
        don=imfilter(don,smoothtype);
        acc=imfilter(acc,smoothtype);
        rfp=imfilter(rfp,smoothtype);
        
        don=don-avgbk;
        doncell=don;
        [m,n]=size(doncell);        
        for j=1:m
            for k=1:n
                if  doncell(j,k)<0
                    doncell(j,k)=0;
                end
            end
        end     
        doncell=doncell./gfpgrad;
        %doncell=circshift(doncell,[rowshift, colshift]); 
        %doncell=imfilter(doncell,smoothtype);
        
        acc=acc-avgbk;
        acccell=acc;        
        for j=1:m
            for k=1:n
                if acccell(j,k)<0
                    acccell(j,k)=0;
                end
            end
        end        
        acccell=acccell./fretgrad;
        %acccell=circshift(acccell,[rowshift, colshift]);                               
        %acccell=imfilter(acccell,smoothtype);
        
        rfp=rfp-avgbk;
        rfpcell=rfp;        
        for j=1:m
            for k=1:n
                if rfpcell(j,k)<0
                    rfpcell(j,k)=0;
                end
            end
        end        
        rfpcell=rfpcell./rfpgrad;
        rfpcell=circshift(rfpcell,[rowshift, colshift]);                       
        %rfpcell=imfilter(rfpcell,smoothtype);
        
        acccellfret=acccell-donleak*doncell-acccross*rfpcell;
        
        for j=1:m
            for k=1:n
                if acccellfret(j,k)<0.0
                    acccellfret(j,k)=0.0;
                end
            end
        end
        
        fret=acccellfret./(rfpcell);

        for j=1:m
            for k=1:n
                %fret(j,k)>=.75 normally but for rebuttal 0.9
                if  fret(j,k)>=.75 || fret(j,k)==inf || fret(j,k)== NaN || doncell(j,k)<intth || rfpcell(j,k)<intth
                    fret(j,k)=0;
                end
            end
        end
        
        %fretraw=fret;
        %fret=imfilter(fret,smoothtype);
        
        if checkintth2==1
           for j=1:m
             for k=1:n
                if  doncell(j,k)<intth2 || rfpcell(j,k)<intth2 || doncell(j,k)>intthmax || rfpcell(j,k)>intthmax
                    fret(j,k)=0.0;
                end
             end
           end
        end           
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(fret, [fretscalemin fretscalemax]);
        colormap(cmap);
        %colorbar;
        axis equal;
        axis off;
        
        fretimg=strcat(opfilepre,'totalcellfret');
        fretimg=strcat(time,'-',fretimg,'.tif');
        %saveas(gca,fretimg,'tiff');                
        tmp=getframe;
        image(tmp.cdata);
        imwrite(tmp.cdata,fretimg);
        
        %%%%%%%%%%%%%% IF FA MASK IS IN THE FOLDER; ONLY FOR REGION FRET VALUES%%%%%%%%%%%%%
        if flag==1
            %figure();
            %[test,rect]=imcrop(don, [50 max(max(don))]); %select region
            %rect=[0,0,imgsizex,imgsizey];
            %fret=imcrop(fret,rect);
            mask=zeros(imgsizex,imgsizey);
            roirect=zeros(roino,4);
            check=1;
            count=0;
            while check==1
                count=count+1;
                %[b rect]=imcrop(don, [50 max(max(don))]); %select region;
                [b rect]=imcrop(don, [50 1000]); %select region;
                rect=uint16(rect);
                
                roirect(count,1)=rect(1);
                roirect(count,2)=rect(2);
                roirect(count,3)=rect(1)+rect(3);
                roirect(count,4)=rect(2)+rect(4);
                
                prompt = {'Select more region if Y enter 1 else 0:'};
                dlg_title = 'More ROI';
                num_lines = 1;
                def = {'1'};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
                if answer{1}~='1'
                    check=0;
                end
            end
            
            for mm=1:count
                for ii=roirect(mm,1):roirect(mm,3)
                    for jj=roirect(mm,2):roirect(mm,4)
                        mask(jj,ii)=1;
                    end
                end
            end
            
            fret=fret.*mask;
        end
            
       if flag==0                               
            %maskfilename=strcat(opfilepre,'Mask of FA-',donfname);   %donfname
            maskfilename='mask-T00001C03Z001.tif';
            %maskfilename=strcat(opfilepre,'mask-',rfpfname);   %donfname
            %maskfilename=strcat(img_direc_base,'FAmask\',donfname);
            mask=imread(maskfilename);
            mask=double(mask);            
            mask=mask./255;
            mask=abs(mask-1.0);            
            fret=fret.*mask;
       end
       if flag==2
           rect=[0,0,imgsizex,imgsizey];
           fret=imcrop(fret,rect);
       end  
       
       if flag==5
           maskfilename=strcat(opfilepre,'Mask of FA-',donfname);
            %maskfilename='Mask of T00001C03Z001.tif';
            %maskfilename=strcat(img_direc_base,'FAmask\',donfname);
            mask=imread(maskfilename);
            mask=double(mask);            
            mask=mask./255;
            %mask=abs(mask-1.0);            
            fret=fret.*mask;
            
           
       end
        
        rect=[0,0,imgsizex,imgsizey];
        
        fretvalfile=strcat(time,'-',opfilepre,'fretvalimage',nonFafilepre,'.txt');
        dlmwrite(fretvalfile,fret,'delimiter','\t','precision',6);
                
        fretlin=fret(:);
        ipos=find(fretlin>0);
        fretlin=fretlin(ipos);
         
        doncelllin=doncell(:);
        %i=find(doncelllin>0);
        doncelllin=doncelllin(ipos);
         
        rfpcelllin=rfpcell(:);
        %i=find(doncelllin>0);
        rfpcelllin=rfpcelllin(ipos);        
        
        if cy5flag==1
            cy5fname= strcat(time,'C03Z001.tif');
            cy5=imread(cy5fname);
            cy5=double(cy5);
            cy5=imfilter(cy5,smoothtype);
            cy5=cy5-avgbk;
            cy5=cy5./cy5grad;
            cy5=circshift(cy5,[rowshift, colshift]);
            cy5lin=cy5(:);
            cy5lin=cy5lin(ipos);
            cy5linfile=strcat(time,'-','cy5int','.txt');
            dlmwrite(cy5linfile,cy5lin,'delimiter','\t','precision',6);
        end
                 
        
        
        fretvalfile=strcat(time,'-',opfilepre,'fretval',nonFafilepre,'.txt');
        dlmwrite(fretvalfile,fretlin,'delimiter','\t','precision',6);
        
        fretvalmeanfile=strcat(time,'-',opfilepre,'fretvalmean',nonFafilepre,'.txt');
        dlmwrite(fretvalmeanfile,mean(fretlin),'delimiter','\t','precision',6);
        
        [bincount,edges]= histcounts(fretlin,edges);
        bincount=bincount';
        bincount=bincount./sum(bincount);
        %plot(edgereal,bincount);
        fretvalhistfile=strcat(time,'-',opfilepre,'frethist',nonFafilepre,'.txt');
        dlmwrite(fretvalhistfile,bincount,'delimiter','\t','precision',6);
        

        figure('units','normalized','outerposition',[0 0 1 1]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %fret=imfilter(fret,smoothtype);
        
        imagesc(fret, [fretscalemin fretscalemax]);
        colormap(cmap);
        %colorbar;
        axis equal;
        axis off;
        fretimg='fret';
        fretimg=strcat(time,'-',opfilepre,nonFafilepre,'fret.tif');
        %saveas(gca,fretimg,'tiff');
        tmp=getframe;
        image(tmp.cdata);
        imwrite(tmp.cdata,fretimg);
        
        fretimg2=strcat(time,'-',opfilepre,nonFafilepre,'fret2.png');
        fret2=uint16(1000*fret);
        imwrite(fret2,fretimg2);
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%         fretbi=zeros(m,n);
%         for j=1:m
%             for k=1:n
%                 if (fret(j,k)>0.0) 
%                     if fret(j,k)<0.075
%                         fretbi(j,k)=100;
%                     else
%                         fretbi(j,k)=200;
%                     end
%                 end
%             end
%         end
%         imagesc(fretbi);
%         colormap(cmap);
%         pause(5);        
%         
%         
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure('units','normalized','outerposition',[0 0 1 1]);
%         doncell=imcrop(doncell,rect);
%         imagesc(doncell, [intth max(max(doncell))]);
%         axis equal;
%         axis off;
%         colormap(cmap);
%         colorbar;
%         donimg=strcat(time,'-','donor.tif');
%         %saveas(gca,donimg,'tiff'); 
%         tmp=getframe;
%         image(tmp.cdata);
%         imwrite(tmp.cdata,donimg);
%       
%         figure('units','normalized','outerposition',[0 0 1 1]);
%         rfpcell=imcrop(rfpcell,rect);
%         imagesc(rfpcell, [intth max(max(rfpcell))]);
%         colormap(cmap);
%         colorbar;
%         axis equal;
%         axis off;   
%         rfpimg=strcat(time,'-','rfp.tif');
%         %saveas(gca,rfpimg,'tiff');
%         tmp=getframe;
%         image(tmp.cdata);
%         imwrite(tmp.cdata,rfpimg);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
         
        %fretlin=(fretlin-min(fretlin))./(max(fretlin)-min(fretlin));
        %doncelllin=(doncelllin-min(doncelllin))./(max(doncelllin)-min(doncelllin));
        %rfpcelllin=(rfpcelllin-min(rfpcelllin))./(max(rfpcelllin)-min(rfpcelllin));
         
        %figure();
        %scatter(rfpcelllin,fretlin);
        %hold on
        %f0=fit(rfpcelllin,fretlin,fittype,'lowess');
        %scatter(rfpcelllin,f0);

%        fretvallinfile=strcat(time,'-','fretvallin','.txt');
%        dlmwrite(fretvallinfile,fretlin,'delimiter','\t','precision',6);
         
        donvallinfile=strcat(time,'-',opfilepre,'donvallin',nonFafilepre,'.txt');
        dlmwrite(donvallinfile,doncelllin,'delimiter','\t','precision',6);
         
        rfpvallinfile=strcat(time,'-',opfilepre,'rfpvallin',nonFafilepre,'.txt');
        dlmwrite(rfpvallinfile,rfpcelllin,'delimiter','\t','precision',6);
    end    
    subdirec1=[subdirec,' ',sprintf('%0d',sub+1)]; 
    pause(2);
end
close all;
beep

        

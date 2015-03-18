%This programcalculates muscle/tendon displacements and velocities
%% **Only moves in horizontal direction
%By Leela Goel

close all
clear
clc

%Input video of interest
%input image name
img_name = input('Name of file? ', 's'); %include file tree and/or extension

%store image
img = VideoReader(img_name);

%specify total number of frames in video
a = get(img, 'NumberOfFrames');
%range of frames to be played
A = [1 a];
vidFrames = read(img, A);

%store frame rate of original image
frame_rate = get(img, 'FrameRate');

%%
for i = 1 : a
    %NEED .CDATA AND .COLORMAP BOTH WHEN CREATING VIDEOS   
    %specify actual image being analyzed
    img_acc(i).cdata = rgb2gray(vidFrames(:,:,:,i)); 
    img_acc(i).colormap = gray;
    
    
    %Crops image
        if i == 1
            fprintf('(1)MEASURE SCALE FIRST, export to workspace and name it pix_per_cm; (2)DONT FORGET TO COPY YOUR IMAGE POSITION!!!!') %need to copy this position by right clicking on selected croped area before double clicking
            pause
    
            %select box for cropping vid
            img_acc2(i).cdata = imcrop(img_acc(i).cdata);
            
            %Determine conversion factor from centimeters to pixels
            %Be sure to measure 1cm from scale on ultrasound vid
            meas_pix_per_cm = imdistline; %measures scale
            pause
            
            % **** EXPORT TO WORKSPACE, CREATE NEW VARIALE 'pix_per_cm' TO
            % USE FOR UNIT CONVERSION
            
            %continue cropping
            im_size = input('Paste your matrix size here'); %hit ctr v or mac equivalent
            img_acc2(i).colormap = gray; 
        else
            img_acc2(i).cdata = imcrop(img_acc(i).cdata, im_size-1);
            img_acc2(i).colormap = gray;
        end
end

%Create Histograms
%{
figure
imhist(img_acc(50).cdata)

figure
imhist(img_acc2(50).cdata)
%}

%%
%delete images and save space
clear img
clear img_acc
clear vidFrames
clear meas_pix_per_cm


%%
%Specify block of interest

frame1 = imshow(img_acc2(1).cdata); %specify first frame for initial region of interest
boi = getrect; %defines coordinates for original block
init_boi = boi; %stores initial boi, for testing purposes

%%
% Determine neighborhood of pixels
pix_per_cm = distance; %gets distance value from exporting line variable
max_pix_disp = round(10 * (1/frame_rate) * pix_per_cm); % 10cm/s is from maximum possible tendon displacement

%%
%define blocks for comparison on next frame 

for i = 1:a-1
    
    %find maximum correlation coefficient
    %then determine new reference box
    %repeat box of interest/neighborhood calculations
  
   
    
    k = 1;
    while k<max_pix_disp
        
        %move box left horizontally
        if (boi(2)) < (im_size(4)) && (boi(1) - (k-1)) < (im_size(3)) %check if box is inbounds of image
            
            ref_boi(1).cdata = imcrop(img_acc2(i).cdata, boi); %define reference image

            c_boi(k).cdata = imcrop(img_acc2(i+1).cdata, [boi(1)- (k-1), boi(2), boi(3), boi(4)]); %get reference block of next frame, fix vertical position from previous motion & block size
            %if j == max_pix_disp
           %corr_coeff(j) = corr_coeff(j-1);
            
           if size(c_boi(k).cdata) == size(ref_boi(1).cdata)
            corr_coeff(k) = corr2(ref_boi(1).cdata, c_boi(k).cdata); %keeps indices from previous matrix of coefficient correlations
           end
        else
            c_boi(k).cdata = c_boi(numel(c_boi)-1).cdata;
            corr_coeff(k) = 0;
        end
        
        if k>1
             if corr_coeff(numel(corr_coeff)) < corr_coeff(numel(corr_coeff) -1)
                
                k_count_cc = k - 1;
                k_ref = k - 1;
                k = max_pix_disp;
                
             end
        end
                 k = k +1;
                 
    end
        
        pix_disp(i) = k_ref; %calculates pixel displacemtn
        boi = [boi(1)- k_ref, boi(2), boi(3), boi(4)]; %creates new corresponding box of interest
end

   
%%
%Convert displacements into centimeters

mus_disp = pix_disp / pix_per_cm; 

%Calculates total displacement

tot_disp = sum(mus_disp)

%%
%Converts displacements into velocities (cm/s)

mus_vel = mus_disp / frame_rate;


%This program will identify and measure muscle fasicle lengths

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

frame1 = imshow(img_acc2(1).cdata) %specify first frame for initial region of interest
frame1
boi = getrect; %defines coordinates for original block

%%
% Determine neighborhood of pixels

max_pix_disp = 10 * (1/frame_rate) * pix_per_cm; % 10cm/s is from maximum possible tendon displacement

%{

%find correlation coefficient between images

for i = 1:a-1
    corr_coeff(i) = corr2(img_acc2(i).cdata, img_acc2(i+1).cdata);
end

%plays cropped vid
%implay(img_acc2)
%}
%%
%convert to grayscale
%img_grsc.cdata = mat2gray(img_acc2.cdata);
%img_grsc.colormap = mat2gray(img_acc2.colormap);


%%
%{

%converts cropped video to double
%use for image intensity array (Strintzis 1997)
for i = 1:a
    img_acc3(i).cdata = double(img_acc2(i).cdata);
    img_acc3(i).colormap = gray;
end

clear img_acc2

%Create Histograms
figure
hist(img_acc3(50).cdata)

%}
%%
%{


%subtract similarities to make only noise and differences shown
%change in motion
%should represent motion field/vector field (Strintzis 1997)
for i = 1:a
    if i == 1
        %Uses first frame from doubled image
        sub_img(i).cdata = img_acc2(i).cdata;
        sub_img(i).colormap = gray;

    else
        %subtracts data from each frame
        sub_img(i).cdata = minus((img_acc2(i).cdata),(img_acc2(i-1).cdata));
        sub_img(i).colormap = gray;

    end
end

clear img_acc3
%Create Histograms
figure
imhist(sub_img(50).cdata)

%}
%%
%{
%add Rayleigh noise to image

for i = 1:a
[M,N] = size(sub_img(i).cdata);

    if i == a
        ray_noi_img(i).cdata = sub_img(i).cdata; %does not add noise to last image
        ray_noi_img(i).colormap = gray;
    else
        ray_noi = imnoise2('rayleigh', M, N, 0, 1);
        ray_noi_img(i).cdata = ray_noi + sub_img(i).cdata; %add noise to image
        ray_noi_img(i).colormap = gray;
    end
end

implay(ray_noi_img)
%Create Histograms
figure
imhist(ray_noi_img(50).cdata)

%}
%%
%{
%determine conditional probabilities 

for i = 1:a-1
    cprob_img(i) = pcond2(ray_noi_img(i).cdata, sub_img(i+1).cdata);
end


%}


%{
%plays subtracted vid
%implay(sub_img);

%applies rayleigh distribution to image
for i = 1:a
    filt_img(i).cdata = adapthisteq(sub_img(i).cdata,'Distribution', 'rayleigh');
    filt_img(i).colormap = [];
end

%clear sub_img

%shows histogram for filtered image
hist(filt_img(50).cdata)
%}    




%%


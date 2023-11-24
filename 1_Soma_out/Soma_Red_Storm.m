%% Get soma filter if necessary. 
%Brain tissue labeling usually comes with unspecific labeling in the soma
%region. This code will identify the soma region based on conventional images of
%selected channel. 

%Commented 1.10.2023. 
clear;clc;

%% Initial setting
exp_folder = 'Y:\Chenghang\ET33_Tigre\20230504_1\analysis\';

%Do you want to apply soma filter?
Soma_filter_use = true;
%Do you want to confine the image using CTB covnentional images? This is
%designed for OPN4 data to confine the analysis in regions that have CTB
%signals. It requires manual processing of CTB conventional images using
%Fourier transformation. 
CTB_Conv_confine = false; 

%Box_thick is used to calculate the boundary of the whole image block. 
%After elastic alginment sometime the cornor is empty. This will help
%including these empty regions in the filter so they won't be calculated in
%the denstiy measurement. 
%Dis_thick is used to expand the identified soma area, , as well as blurring of the bouder.
%Channel: R=1,G=2,B=3.
%gauss_id: the gaussian blur parameter. Higher if the soma staining is not
%prominent.
Box_thick = 6;
Dis_thick = 40;
channel = 1;
gauss_id = 3;

path = [exp_folder  'elastic_align_rescale/'];
stormpath = [path 'storm_merged/'];
stormfiles = [dir([stormpath '*.tif']) dir([stormpath '*.png'])];
num_images = numel(stormfiles);
infos = imfinfo([stormpath stormfiles(1,1).name]);
convpath = [path 'conv_merged\'];
stormpath = [path 'storm_merged\'];
convfiles = [dir([convpath '*.tif']) dir([convpath '*.png'])];

%convpath_thre is used when using Fouriour method manually
convpath_thre = [path 'conv_488_thre_1\'];
conv_thre_files = [dir([convpath_thre '*.tif']) dir([convpath_thre '*.png'])];

mkdir([exp_folder 'Result']);
mkdir([exp_folder 'Result\1_soma']);
mkdir([exp_folder 'Result\0_BD'])
outpath = [exp_folder 'Result\1_Soma\'];
voxels=[(15.5), (15.5), 70];

BYs = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
BD = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');

% Find boundary of each image section
parfor k = 1:num_images
    disp(k)
    A = imread([convpath convfiles(k,1).name]);
    C = A(:,:,channel) * 255;
    BD(:,:,k) = C;
end
BD = logical(BD);
%
%% Find soma regions
%You can use convfiles/stormfiles to identify the soma regions. 
%Manually change the t_use value may be necessary. 
if Soma_filter_use == true
    parfor k = 1:num_images
        A = imread([convpath convfiles(k,1).name]);
        %A = imread([stormpath convfiles(k,1).name]);
        BYs(:,:,k) = A(:,:,channel);
    end

    lo_int = 0;
    hi_int = stretchlim(BYs(:),0.003);
    hi_int = hi_int(2);
    for i = 1:num_images
        BYs(:,:,i) = imadjust(BYs(:,:,i),[lo_int,hi_int]);
    end

    gausspix = (gauss_id);
    bg2 = zeros(size(BYs),'uint8');
    parfor k=1:num_images
        disp(k)
        bg2(:,:,k) = (imfilter(BYs(:,:,k),fspecial('gaussian',gausspix,gausspix),'same','replicate'));
    end
    %
    %
    %
    %Change this value to test the best soma filter. 
    t_use = 0.55;
    
    disp('making CG')
    CG = false(size(bg2));
    parfor k=1:size(bg2,3)
        CG(:,:,k) = im2bw(bg2(:,:,k), t_use);
    end
    %
    disp('making CCG')
    CCG = bwconncomp(CG,26);
    %clear CG
    disp('making statsG')

    statsG = regionprops(CCG,BYs,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');
    statsGgauss = regionprops(CCG,bg2,'Area','PixelIdxList','PixelValues','PixelList','WeightedCentroid');

    %
    new_G = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
    %
    % the size threshold here might need modification sometimes.
    statsG_temp = statsGgauss(find(log(1+[statsGgauss.Area])> 12)); %11.5
    numel(statsG_temp)
    %
    for i = 1:size(statsG_temp,1)
        % for i = 1:10
        disp(i)
        PixelList = statsG_temp(i).PixelList;
        PixelValues = statsG_temp(i).PixelValues;
        for j = 1:size(PixelList,1)
            y = PixelList(j,1);
            x = PixelList(j,2);
            z = PixelList(j,3);
            new_G(x,y,z) = 1;
        end
    end
    %
    %Write the images for a quick check. 
    imwrite(new_G(:,:,1) * 255,[outpath 'Soma_test_1_' char(string(t_use)) '.tif']);
    imwrite(new_G(:,:,20) * 255,[outpath 'Soma_test_20_' char(string(t_use)) '.tif']);
else
    %WHen not use the soma filter, get an empty new_G. 
    new_G = zeros(infos.Height, infos.Width, num_images,'uint8');
end
%% Expand and fill the soma region. 
new_G = logical(new_G);
soma_filter = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
parfor i = 1:num_images
    Dg = bwdist(new_G(:,:,i));
    curr = Dg<=(Dis_thick);
    curr = imfill(curr,4,'holes');
    soma_filter(:,:,i) = curr;
end
%% Confine the region using Conventional CTB signals. 
% Run these code step by step to determin the best parameters. 
if CTB_Conv_confine == true
    %add the CTB-dense area filter here: Get the centroid list
    CTB_dense_Area = zeros(infos.Height, infos.Width, num_images);
    CTB_dense_Area = logical(CTB_dense_Area);
    parfor i = 1:num_images
        temp = imread([convpath_thre sprintf('%03d',i) '.png']);
        temp = logical(temp);
        CTB_dense_Area(:,:,i) = temp;
    end
    CTB_conn = bwconncomp(CTB_dense_Area,26);
    CTB_conn_props = regionprops(CTB_conn,CTB_dense_Area,'PixelIdxList','PixelValues','WeightedCentroid');
    Centroids = zeros(numel(CTB_conn_props),3);
    for i = 1:numel(CTB_conn_props)
        temp = CTB_conn_props(i).WeightedCentroid;
        Centroids(i,1) = temp(1);
        Centroids(i,2) = temp(2);
        Centroids(i,3) = temp(3);
    end
    %
    %Find the distance from the centered point to all other points.
    Centroids(:,1) = Centroids(:,1).*0.0155;
    Centroids(:,2) = Centroids(:,2).*0.0155;
    Centroids(:,3) = Centroids(:,3).*0.0700;
    Tot_intensity = zeros(numel(CTB_conn_props),1);
    for i = 1:numel(CTB_conn_props)
        Tot_intensity(i) = sum(CTB_conn_props(i).PixelValues);
    end
    Mean_centroid(1) = sum(Centroids(:,1).*Tot_intensity) / sum(Tot_intensity);
    Mean_centroid(2) = sum(Centroids(:,2).*Tot_intensity) / sum(Tot_intensity);
    Mean_centroid(3) = sum(Centroids(:,3).*Tot_intensity) / sum(Tot_intensity);

    dist_matrix = sqrt((Centroids(:,1) - Mean_centroid(1)).^2 + (Centroids(:,2) - Mean_centroid(2)).^2 + (Centroids(:,3) - Mean_centroid(3)).^2);
    figure;plot(sort(dist_matrix));
    %%
    %Manually select a threhsold.
    Dist_thre = 60;
    Centroids_within = Centroids(dist_matrix<=Dist_thre,:);
    figure;
    subplot(1,2,1);
    plot(Centroids(:,1),Centroids(:,2),'.');
    xlim([0,140]);ylim([0,140])
    subplot(1,2,2);
    plot(Centroids(:,1),Centroids(:,2),'.');hold on;
    plot(Centroids_within(:,1),Centroids_within(:,2),'r.');
    xlim([0,140]);ylim([0,140])
    %%
    %Find the convex hull and apply it to the filter (CTB_Conv_convex).
    CTB_conn_props_within = CTB_conn_props(dist_matrix<=Dist_thre,:);
    CTB_dense_Area = zeros(infos.Height, infos.Width, num_images);
    for i = 1:numel(CTB_conn_props_within)
        PixelValues = CTB_conn_props_within(i).PixelValues;
        PixelIdxList = CTB_conn_props_within(i).PixelIdxList;
        CTB_dense_Area(PixelIdxList) = PixelValues;
    end
    for i = 1:num_images
        CTB_Conv_convex(:,:,i) = bwconvhull(logical(CTB_dense_Area(:,:,i)));
    end
end
%% Finalize the filter. 
%Get neuropil regions and all selected regions. 
if CTB_Conv_confine == true
    CTB_Conv_convex = uint8(CTB_Conv_convex);
    neuropil_filter = ((logical(CTB_Conv_convex) - logical(soma_filter)) > 0);
else
    CTB_Conv_convex = BD;
    neuropil_filter = ((CTB_Conv_convex - logical(soma_filter)) > 0);
end

Area_used = zeros(num_images,1);
for i = 1:num_images
    Area_used(i) = numel(find(CTB_Conv_convex(:,:,i)));
end

Area_neuropil = zeros(num_images,1);
for i = 1:num_images
    Area_neuropil(i) = numel(find(neuropil_filter(:,:,i)));
end

for i = 1:num_images
    %Write the neuropil region. 
    imwrite(double(neuropil_filter(:,:,i)),[outpath 'F_' sprintf('%03d',i) '.tif']);
    %Write all the regions. 
    imwrite(double(CTB_Conv_convex(:,:,i)),[exp_folder 'Result/0_BD/' 'A2_' sprintf('%03d',i) '.tif']);
end
%% Quantification save out
%Save the voxel number to a seperate folder. 
save_path = 'Y:\Chenghang\ET33_Tigre\20230504_1\analysis\Result\Quantification\';
temp = Area_used;
save([save_path 'Area_used.txt'],'temp','-ascii','-double');
temp = Area_neuropil;
save([save_path 'Area_neuropil.txt'],'temp','-ascii','-double');
save([exp_folder 'Result\1_soma\Soma_stats.mat'],'statsG_temp','-v7.3');
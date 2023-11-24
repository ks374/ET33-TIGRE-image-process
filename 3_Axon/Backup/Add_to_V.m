%% add_to_R.m
%Code explanation see 'Add_to_G.m'
%%
clear;clc
%
base_path = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\';
exp_folder = [base_path 'analysis\'];
path = [exp_folder  'elastic_align\'];
mergedpath = [path 'storm_merged/'];
%
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);
disp(num_images)
%set voxel size in nm
voxel=[15.5, 15.5, 70];

outpath = [exp_folder '\Result\3_Vglut2\'];
outpath_syn = [exp_folder '\Result\'];
%
load([outpath_syn 'R_paired_2.mat']);
statsRwater = statsRwater_sss;
clear statsGwater
%
% load([outpath_syn 'G_paired_2.mat']);
% statsGwater = statsGwater_sss;
%
statsSwater = statsRwater;
%
load([outpath 'statsV2sw10.mat']);
statsVwater = statsGa2s;

%
clear BP BG bg2
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
BP2 = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsVwater)
    disp(int2str(i))
    BP(statsVwater(i).PixelIdxList)=statsVwater(i).PixelValues;
end
%
%BP2 = BP;
for i = 1:numel(statsSwater)
    disp(int2str(i))
    BP2(statsSwater(i).PixelIdxList)=statsSwater(i).PixelValues;
end
%
disp('loading data')
%
%if not present, load in stats list

%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsVwater):-1:1
    disp(i)
    
    statsVwater(i).tints_p140 = [];
    statsVwater(i).volume_p140 = [];
    statsVwater(i).area_p140 = [];
    statsVwater(i).WeightedCentroid_p140 = [];
    statsVwater(i).statxtureP_140all = [];
    statsVwater(i).statxtureP_140pos = [];
    centG_p140(i,1:3) = 0;
    tintsG_p140(i) = 0;
    volumeG_p140(i) =0;
    
    statsVwater(i).statxture_all = [];
    statsVwater(i).statxture_pos = [];
    statxture_all(i,1:6) = 0;
    statxture_pos(i,1:6) = 0;
    
    %try preproccessing the BP varialbe in a for loop
    minpix = min(statsVwater(i).PixelList);  maxpix = max(statsVwater(i).PixelList);
    min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
    max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
    if min1 < 1; min1=1; end
    if min2 < 1; min2=1; end
    if min3 < 1; min3=1; end
    if max1 > info.Width; max1=info.Width; end
    if max2 > info.Height; max2=info.Height; end
    if max3 > num_images; max3=num_images; end
    BPp(i).mat = BP(min2:max2,min1:max1,min3:max3);
    BP2p(i).mat = BP2(min2:max2,min1:max1,min3:max3);
end
%

parfor jj=1:numel(statsVwater)
    %
    disp(jj)
    minpix = min(statsVwater(jj).PixelList);  maxpix = max(statsVwater(jj).PixelList);
    min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
    max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
    if min1 < 1; min1=1; end
    if min2 < 1; min2=1; end
    if min3 < 1; min3=1; end
    if max1 > info.Width; max1=info.Width; end
    if max2 > info.Height; max2=info.Height; end
    if max3 > num_images; max3=num_images; end
    
    size1 = max1-min1 + 1; size2 = max2-min2 + 1; size3 = max3-min3 + 1;
    curr2 = false(size2, size1, size3);
    
    for j=1: numel(statsVwater(jj).PixelList(:,1))
        curr2(statsVwater(jj).PixelList(j,2)-min2+1, ...
            statsVwater(jj).PixelList(j,1)-min1+1, ...
            statsVwater(jj).PixelList(j,3)-min3+1)= 1;
    end
    curr1a = BPp(jj).mat;
    curr1b = BP2p(jj).mat;
    Dg = bwdistsc1(curr2,[voxel(1),voxel(2),voxel(3)]);
    %
    
    size2= 140;
    curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
    [y,x,z] = ind2sub(size(curr4b),find(curr4b));
    PixelList_m140 = [x,y,z];
    RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
    for kk=1:numel(PixelList_m140(:,1))
        RPixelValues_m140(kk) = curr1b(PixelList_m140(kk,2),...
            PixelList_m140(kk,1), PixelList_m140(kk,3));
    end
    
    statsVwater(jj).tints_p140 = sum([RPixelValues_m140]);
    statsVwater(jj).area_p140 = numel([RPixelValues_m140]);
    statsVwater(jj).volume_p140 = numel([RPixelValues_m140]>0);
    statsVwater(jj).WeightedCentroid_p140(1) = sum([PixelList_m140(:,1)].*...
        double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
    statsVwater(jj).WeightedCentroid_p140(2) = sum([PixelList_m140(:,2)].*...
        double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
    statsVwater(jj).WeightedCentroid_p140(3) =  sum([PixelList_m140(:,3)].*...
        double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
    %also add these statxture,
    %    statsVwater(jj).statxtureP_140all = statxture([RPixelValues_m140]);
    %    posvals2 = [RPixelValues_m140]>0;
    %    statsVwater(jj).statxtureP_70pos = statxture([RPixelValues_m140(posvals2)]);
    centG_p140(jj,:) = statsVwater(jj).WeightedCentroid_p140;
    tintsG_p140(jj) = statsVwater(jj).tints_p140;
    volumeG_p140(jj) =statsVwater(jj).volume_p140;
    
    
    
    
    %    statsVwater(jj).statxture_all = statxture([statsVwater(jj).PixelValues]);
    %    posvals = [statsVwater(jj).PixelValues]>0;
    %    statsVwater(jj).statxture_pos = statxture([statsVwater(jj).PixelValues(posvals)]);
    %    statxture_all(jj,:) = statsVwater(jj).statxture_all;
    %    statxture_pos(jj,:) = statsVwater(jj).statxture_pos;
end
%
save([outpath 'add_to_statsVw10_edges.mat'],'centG_p140','tintsG_p140',...)
    'volumeG_p140','statxture_all','statxture_pos','-v7.3')
save([outpath 'statsV2w10_edges_plus.mat'],'statsVwater','-v7.3')
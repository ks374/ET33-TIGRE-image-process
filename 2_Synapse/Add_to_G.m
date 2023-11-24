%% add_to_statsG
%Based on identifiied G and R clusters, for each G cluster, find signal of
%R cluster within a certain range. The range can be set in 'size2' in the
%unit of nm. Default 140. 
%%
clear;clc;
%
base_path = 'Y:\Chenghang\ET33_Tigre\20230504_1\';
exp_folder = [base_path 'analysis\'];
path = [exp_folder  'elastic_align_rescale\'];
mergedpath = [path 'storm_merged\'];
%
mergedfiles = [dir([mergedpath '*.png']) dir([mergedpath '*.tif'])];
num_images = numel(mergedfiles);
info = imfinfo([mergedpath mergedfiles(1,1).name]);
disp(num_images)
%set voxel size in nm
voxel=[15.5, 15.5, 70];

outpath = [exp_folder 'Result\'];
%
load([outpath 'statsR2sw10.mat']);
statsRwater = statsGa2s;
clear statsGwater
%
clear statsGa2s statsGwater
load([outpath 'statsG2sw10.mat']);
statsGwater = statsGa2s;
%
% num_images = num_images - 2;
clear BP BG bg2
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
BP2 = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsGwater)
    disp(int2str(i))
    BP(statsGwater(i).PixelIdxList)=statsGwater(i).PixelValues;  
end
%
%BP2 = BP;
for i = 1:numel(statsRwater)
    disp(int2str(i))
    BP2(statsRwater(i).PixelIdxList)=statsRwater(i).PixelValues;
end
%
disp('loading data')
%
%if not present, load in stats list

%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsGwater):-1:1
    disp(i)
    

statsGwater(i).tints_p140 = [];
statsGwater(i).volume_p140 = [];
statsGwater(i).area_p140 = [];
statsGwater(i).WeightedCentroid_p140 = [];
statsGwater(i).statxtureP_140all = [];
statsGwater(i).statxtureP_140pos = [];
centG_p140(i,1:3) = 0;
tintsG_p140(i) = 0;
volumeG_p140(i) =0;

statsGwater(i).statxture_all = [];
statsGwater(i).statxture_pos = [];
statxture_all(i,1:6) = 0;
statxture_pos(i,1:6) = 0;
%try preproccessing the BP varialbe in a for loop
   minpix = min(statsGwater(i).PixelList);  maxpix = max(statsGwater(i).PixelList);
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

parfor jj=1:numel(statsGwater)
   %
   disp(jj)
   minpix = min(statsGwater(jj).PixelList);  maxpix = max(statsGwater(jj).PixelList);
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

   for j=1: numel(statsGwater(jj).PixelList(:,1))
       curr2(statsGwater(jj).PixelList(j,2)-min2+1, ...
          statsGwater(jj).PixelList(j,1)-min1+1, ...
          statsGwater(jj).PixelList(j,3)-min3+1)= 1;  
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

   statsGwater(jj).tints_p140 = sum([RPixelValues_m140]);
   statsGwater(jj).area_p140 = numel([RPixelValues_m140]);
   statsGwater(jj).volume_p140 = numel([RPixelValues_m140]>0);
   statsGwater(jj).WeightedCentroid_p140(1) = sum([PixelList_m140(:,1)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsGwater(jj).WeightedCentroid_p140(2) = sum([PixelList_m140(:,2)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsGwater(jj).WeightedCentroid_p140(3) =  sum([PixelList_m140(:,3)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
%also add these statxture, 
%    statsGwater(jj).statxtureP_140all = statxture([RPixelValues_m140]);
%    posvals2 = [RPixelValues_m140]>0;
%    statsGwater(jj).statxtureP_70pos = statxture([RPixelValues_m140(posvals2)]);
   centG_p140(jj,:) = statsGwater(jj).WeightedCentroid_p140;
   tintsG_p140(jj) = statsGwater(jj).tints_p140;
   volumeG_p140(jj) =statsGwater(jj).volume_p140;
   
   
 
   
%    statsGwater(jj).statxture_all = statxture([statsGwater(jj).PixelValues]);
%    posvals = [statsGwater(jj).PixelValues]>0;
%    statsGwater(jj).statxture_pos = statxture([statsGwater(jj).PixelValues(posvals)]);
%    statxture_all(jj,:) = statsGwater(jj).statxture_all;  
%    statxture_pos(jj,:) = statsGwater(jj).statxture_pos;
end
%
save([outpath 'add_to_statsGw10_edges.mat'],'centG_p140','tintsG_p140',...)
    'volumeG_p140','statxture_all','statxture_pos','-v7.3')
save([outpath 'statsG2w10_edges_plus.mat'],'statsGwater','-v7.3')
%% Add_to.m
% From B clusters to G clusters. 
%%
clear;clc
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

Syn_path = [base_path '\analysis\Result\'];
V_path = [base_path '\analysis\Result\3_Vglut2\'];
mkdir([base_path '\analysis\Result\5_V_Syn\']);
outpath = [base_path '\analysis\Result\5_V_Syn\'];
%
load([Syn_path 'R_paired_2.mat']);
statsRwater = statsRwater_sss;

load([V_path 'V_paired.mat']);
statsVwater = statsVwater_ss;

%
clear BP BG bg
disp('allocating arrays')
BP = zeros(info.Height, info.Width, num_images,'uint8');
BP2 = zeros(info.Height, info.Width, num_images,'uint8');
for i = 1:numel(statsRwater)
    disp(int2str(i))
    BP(statsRwater(i).PixelIdxList)=statsRwater(i).PixelValues;  
end
%
%BP2 = BP;
for i = 1:numel(statsVwater)
    disp(int2str(i))
    BP2(statsVwater(i).PixelIdxList)=statsVwater(i).PixelValues;  
end
%
disp('loading data')
%
%if not present, load in stats list

%make categories for new variables in parfor loop and slice BP for use in
%parfor loop
for i=numel(statsRwater):-1:1
    disp(i)

statsRwater(i).tints_p140 = [];
statsRwater(i).volume_p140 = [];
statsRwater(i).area_p140 = [];
statsRwater(i).WeightedCentroid_p140 = [];
statsRwater(i).statxtureP_140all = [];
statsRwater(i).statxtureP_140pos = [];
centG_p140(i,1:3) = 0;
tintsG_p140(i) = 0;
volumeG_p140(i) =0;

statsRwater(i).statxture_all = [];
statsRwater(i).statxture_pos = [];
statxture_all(i,1:6) = 0;
statxture_pos(i,1:6) = 0;

%try preproccessing the BP varialbe in a for loop
   minpix = min(statsRwater(i).PixelList);  maxpix = max(statsRwater(i).PixelList);
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

parfor jj=1:numel(statsRwater)
%for jj=500:600
   %
   disp(jj)
   minpix = min(statsRwater(jj).PixelList);  maxpix = max(statsRwater(jj).PixelList);
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

   for j=1: numel(statsRwater(jj).PixelList(:,1))
       curr2(statsRwater(jj).PixelList(j,2)-min2+1, ...
          statsRwater(jj).PixelList(j,1)-min1+1, ...
          statsRwater(jj).PixelList(j,3)-min3+1)= 1;  
   end
   curr1a = BPp(jj).mat;
   curr1b = BP2p(jj).mat;
   Dg = bwdistsc1(curr2,[voxel(1),voxel(2),voxel(3)]);
%

   size2= 70;
   curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
   [y,x,z] = ind2sub(size(curr4b),find(curr4b));
   PixelList_m140 = [x,y,z];
   RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
   for kk=1:numel(PixelList_m140(:,1))
      RPixelValues_m140(kk) = curr1b(PixelList_m140(kk,2),...
          PixelList_m140(kk,1), PixelList_m140(kk,3));
   end

   statsRwater(jj).tints_p140 = sum([RPixelValues_m140]);
   statsRwater(jj).area_p140 = numel([RPixelValues_m140]);
   statsRwater(jj).volume_p140 = numel([RPixelValues_m140]>0);
   statsRwater(jj).WeightedCentroid_p140(1) = sum([PixelList_m140(:,1)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsRwater(jj).WeightedCentroid_p140(2) = sum([PixelList_m140(:,2)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
   statsRwater(jj).WeightedCentroid_p140(3) =  sum([PixelList_m140(:,3)].*...
           double([RPixelValues_m140]))/(sum([RPixelValues_m140]));
%also add these statxture, 
%    statsRwater(jj).statxtureP_140all = statxture([RPixelValues_m140]);
%    posvals2 = [RPixelValues_m140]>0;
%    statsRwater(jj).statxtureP_70pos = statxture([RPixelValues_m140(posvals2)]);
   centG_p140(jj,:) = statsRwater(jj).WeightedCentroid_p140;
   tintsG_p140(jj) = statsRwater(jj).tints_p140;
   volumeG_p140(jj) =statsRwater(jj).volume_p140;
  

   
%    statsRwater(jj).statxture_all = statxture([statsRwater(jj).PixelValues]);
%    posvals = [statsRwater(jj).PixelValues]>0;
%    statsRwater(jj).statxture_pos = statxture([statsRwater(jj).PixelValues(posvals)]);
%    statxture_all(jj,:) = statsRwater(jj).statxture_all;
%    statxture_pos(jj,:) = statsRwater(jj).statxture_pos;
end
%
save([outpath 'add_to_statsRw10_edges_Vglut2.mat'],'centG_p140','tintsG_p140',...)
    'volumeG_p140','statxture_all','statxture_pos','-v7.3');
save([outpath 'statsR2w10_edges_plus_Vglut2.mat'],'statsRwater','-v7.3');
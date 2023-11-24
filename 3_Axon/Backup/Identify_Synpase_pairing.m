%% Pair B clusters to G
%This is a simplified version of synapse pairing. 
%%
base = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\';
base_folder = [base 'analysis\Result\3_Vglut2\'];
outpath = base_folder;
outpath2 = [base 'analysis\Result\'];
%
load([outpath2 'G_paired_2.mat']);
voxel = [15.5 15.5 70];

for i =1:numel(statsGwater_sss)
    volumeRs(i,1) = statsGwater_sss(i).Volume1_0;
end

%
load([base_folder 'statslistV2sw10.mat']);
load([base_folder 'statslistV2nsw10.mat']);

volumeGs = sizeshape_matGa2s(:,19);
%
load([base_folder 'add_to_statsVw10_edges.mat'],'tintsG_p140');
mints_g70s = (([tintsG_p140])./[volumeGs]');
%
pairedg_idx = find(mints_g70s);
disp('The number of selected Vgltu2 cluster: ')
numel(pairedg_idx)/numel(mints_g70s)

V_paired_id = zeros(numel(mints_g70s),1);
V_paired_id(pairedg_idx) = 1;
V_paired_id = logical(V_paired_id);
pairedg_idx = V_paired_id;
save([outpath 'V_paired_id.mat'],'pairedg_idx');
%
load([outpath 'statsV2sw10.mat']);
statsVwater_ss = statsGa2s(pairedg_idx);
statsVwater_sn = statsGa2s(~pairedg_idx);
%% 0.1% saturated image restore.
%From the 0.3% saturated clusters, get corresponded information from 0.1%
%saturated images. 
%Don't need this if 0.3% saturated images were used from the beginning. 
Storm_path = [base 'analysis\elastic_align\storm_merged\'];
files = [dir([Storm_path '*.tif'])];
infos = imfinfo([Storm_path files(1,1).name]);
num_images = numel(files);
Storm_images_B = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');

parfor i = 1:num_images
    Storm_image = imread([Storm_path sprintf('%03d',i) '.tif']);
    Storm_images_B(:,:,i) = Storm_image(:,:,3);
end
statsVwater_ss_01 = [];
for i = 1:numel(statsVwater_ss)
    stats_temp = stats_copy(statsVwater_ss(i),Storm_images_B);
    statsVwater_ss_01 = cat(2,statsVwater_ss_01,stats_temp);
end
statsVwater_sn_01 = [];
for i = 1:numel(statsVwater_sn)
    stats_temp = stats_copy(statsVwater_sn(i),Storm_images_B);
    statsVwater_sn_01 = cat(2,statsVwater_sn_01,stats_temp);
end

clear statsVwater_ss statsVwater_sn
statsVwater_ss = statsVwater_ss_01;
statsVwater_sn = statsVwater_sn_01;
%
save([base_folder 'V_paired.mat'],'statsVwater_ss','statsVwater_sn');
%% Quantificaton. 
save_path = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\analysis\Result\Quantification\';
temp = [statsVwater_ss.Volume1_0];
temp = temp';
save([save_path 'statsVwater_ss.Area.txt'],'temp','-ascii','-double');
temp = [statsVwater_ss.TintsG];
temp = temp';
save([save_path 'statsVwater_ss.TintsG.txt'],'temp','-ascii','-double');

%% G pairing. 
%A simplified code to find R clusters that are paired with G. 
%%
base = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\';
base_folder = [base '\analysis\Result\5_V_Syn\'];
outpath = base_folder;
%
load([base_folder 'statsG2w10_edges_plus_Vglut2.mat']);
voxel = [15.5 15.5 70];
%
for i =1:numel(statsGwater)
    volumeGs(i,1) = statsGwater(i).Volume1_0;
end

%
load([base_folder 'add_to_statsGw10_edges_Vglut2.mat'],'tintsG_p140');
mints_g70s = (([tintsG_p140])./[volumeGs]');
%
pairedg_idx = find(mints_g70s);
disp('The number of filtered out suspected Vgltu2 cluster: ')
disp(numel(pairedg_idx)/numel(mints_g70s));

G_paired_GV_id = zeros(numel(mints_g70s),1);
G_paired_GV_id(pairedg_idx) = 1;
G_paired_GV_id = logical(G_paired_GV_id);
pairedg_idx = G_paired_GV_id;
save([outpath 'G_paired_GV_id.mat'],'pairedg_idx');
%
statsGwater_ssss = statsGwater(pairedg_idx);
statsGwater_sssn = statsGwater(~pairedg_idx);

save([base_folder 'G_paired_3.mat'],'statsGwater_sssn','statsGwater_ssss');

%%
save_path = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\analysis\Result\Quantification\';
temp = [statsGwater_ssss.Volume1_0];
temp = temp';
save([save_path 'statsGwater_ssss.Area.txt'],'temp','-ascii','-double');
temp = [statsGwater_ssss.TintsG];
temp = temp';
save([save_path 'statsGwater_ssss.TintsG.txt'],'temp','-ascii','-double');
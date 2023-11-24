%% G pairing. 
%A simplified code to find G clusters that are paired with B. 
%%
base = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\';
base_folder = [base 'analysis\Result\5_V_Syn\'];
outpath = base_folder;
%
load([base_folder 'statsR2w10_edges_plus_Vglut2.mat']);
voxel = [15.5 15.5 70];
%
for i =1:numel(statsRwater)
    volumeGs(i,1) = statsRwater(i).Volume1_0;
end

%
load([base_folder 'add_to_statsRw10_edges_Vglut2.mat'],'tintsG_p140');
mints_g70s = (([tintsG_p140])./[volumeGs]');
%
pairedg_idx = find(mints_g70s);
disp('The number of filtered out suspected Vgltu2 cluster: ')
numel(pairedg_idx)/numel(mints_g70s)

R_paired_V_id = zeros(numel(mints_g70s),1);
R_paired_V_id(pairedg_idx) = 1;
R_paired_V_id = logical(R_paired_V_id);
pairedg_idx = R_paired_V_id;
save([outpath 'R_paired_V_id.mat'],'pairedg_idx');
%
statsRwater_ssss = statsRwater(pairedg_idx);
statsRwater_sssn = statsRwater(~pairedg_idx);
save([base_folder 'R_paired_3.mat'],'statsRwater_sssn','statsRwater_ssss');
%%
save_path = 'Y:\Chenghang\4_Color\Raw\1.11.2021_B2P8C_B_V2\analysis\Result\Quantification\';
temp = [statsRwater_ssss.Volume1_0];
temp = temp';
save([save_path 'statsRwater_ssss.Area.txt'],'temp','-ascii','-double');
temp = [statsRwater_ssss.TintsG];
temp = temp';
save([save_path 'statsRwater_ssss.TintsG.txt'],'temp','-ascii','-double');

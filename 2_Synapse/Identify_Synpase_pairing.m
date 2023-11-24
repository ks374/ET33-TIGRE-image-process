%% Paire R and G clusters. 
% Require output from Add_to_G and Add_to_R. 
clear;clc
%
exp_folder = 'Y:\Chenghang\ET33_Tigre\20230504_1\';
base_folder = [exp_folder '\analysis\Result\'];
outpath = base_folder;
%
load([base_folder 'statslistR2sw10.mat']);
load([base_folder 'statslistR2nsw10.mat']);
%
centRa2ns = centGa2ns;
centRa2s = centGa2s;
rcentRa2ns = rcentGa2ns;
rcentRa2s = rcentGa2s;
sizeshape_matRa2ns = sizeshape_matGa2ns;
sizeshape_matRa2s = sizeshape_matGa2s;
clear centGa2ns centGa2s rcentGa2ns rcentGa2s sizeshape_matGa2ns sizeshape_matGa2s

load([base_folder 'statslistG2sw10.mat']);
load([base_folder 'statslistG2nsw10.mat']);

%
sizeshape_matGns = sizeshape_matGa2ns;%(rcentGa2ns(:,2)<46000,:);
centGa2ns2 = centGa2ns;%(rcentGa2ns(:,2)<46000,:);
centGns = rcentGa2ns;%(rcentGa2ns(:,2)<46000,:);

sizeshape_matGs = sizeshape_matGa2s;%(rcentGa2s(:,2)<46000,:);
centGa2s2 = centGa2s;%(rcentGa2s(:,2)<46000,:);
centGs = rcentGa2s;%(rcentGa2s(:,2)<46000,:);
%
sizeshape_matRns = sizeshape_matRa2ns;%(rcentPa2ns(:,2)<46000,:);
centRa2ns2 = centRa2ns;%(rcentPa2ns(:,2)<46000,:);
centRns = rcentRa2ns;%(rcentPa2ns(:,2)<46000,:);

sizeshape_matRs = sizeshape_matRa2s;%(rcentPa2s(:,2)<46000,:);
centRa2s2 = centRa2s;%(rcentPa2s(:,2)<46000,:);
centRs = rcentRa2s;%(rcentPa2s(:,2)<46000,:);

centR_all = cat(1, centRs, centRns);
volumeRs = sizeshape_matRs(:,19);
centG_all = cat(1, centGns, centGs);
volumeGs = sizeshape_matGs(:,19);

%
expfolder = [exp_folder '\analysis\elastic_align_rescale\storm_merged\'];
files = [dir([expfolder '*.tif']) dir([expfolder '*.png'])];
infos = imfinfo([expfolder files(1,1).name]);
num_images = numel(files);
size1 = infos.Width;
size2 = infos.Height;
sizeR = size(centRs,1);
sizeRn = size(centRns,1);
sizeRall = size(centR_all,1);
sizeG = size(centGs,1);
sizeGn = size(centGns,1);
sizeGall = size(centG_all,1);
%
voxel = [15.5,15.5,70];
centRs_rand = cat(2,randi(size1,sizeR,1).*voxel(1),randi(size2,sizeR,1).*voxel(2),randi(num_images,sizeR,1).*voxel(3));
centRall_rand = cat(2,randi(size1,sizeRall,1).*voxel(1),randi(size2,sizeRall,1).*voxel(2),randi(num_images,sizeRall,1).*voxel(3));
centGs_rand = cat(2,randi(size1,sizeG,1).*voxel(1),randi(size2,sizeG,1).*voxel(2),randi(num_images,sizeG,1).*voxel(3));
centGall_rand = cat(2,randi(size1,sizeGall,1).*voxel(1),randi(size2,sizeGall,1).*voxel(2),randi(num_images,sizeGall,1).*voxel(3));
%
nn_Gall_Rs = zeros(size(centG_all,1),1);
nn_Grand_Rs = zeros(size(centG_all,1),1);
nn_Gs_Rs = zeros(size(centGs,1),1);
nn_Gns_Rs = zeros(size(centGns,1),1);
%nn_Gall_Pall = zeros(size(centG_all,1),1);
%nn_Grand_Pall = zeros(size(centG_all,1),1);
%nn_Gs_Pall = zeros(size(centGs,1),1);
%nn_Gns_Pall = zeros(size(centGns,1),1);

nn_Rall_Gs = zeros(size(centR_all,1),1);
nn_Rrand_Gs = zeros(size(centR_all,1),1);
nn_Rs_Gs = zeros(size(centRs,1),1);
nn_Rns_Gs = zeros(size(centRns,1),1);
%nn_Pall_Gall = zeros(size(centP_all,1),1);
%nn_Prand_Gall = zeros(size(centP_all,1),1);
%nn_Ps_Gall = zeros(size(centPs,1),1);
%nn_Pns_Gall = zeros(size(centPns,1),1);
%
disp('startG1')
parfor i=1:size(centG_all,1)
nn_Gall_Rs(i) = min(pdist2(centG_all(i,:),centRs));
nn_Grand_Rs(i) = min(pdist2(centGall_rand(i,:),centRs_rand));
%nn_Gall_Pall(i) = min(pdist2(centG_all(i,:),centP_all));
%nn_Grand_Pall(i) = min(pdist2(centGall_rand(i,:),centP_all));
end
%
disp('startG2')
parfor i=1:size(centGs,1)
   nn_Gs_Rs(i) = min(pdist2(centGs(i,:),centRs));
%   nn_Gs_Pall(i) = min(pdist2(centGs(i,:),centP_all));
end
disp('startG3')
parfor i=1:size(centGns,1)
   nn_Gns_Rs(i) = min(pdist2(centGns(i,:),centRs));
%   nn_Gns_Pall(i) = min(pdist2(centGns(i,:),centP_all));
end
%
disp('startP1')
parfor i=1:size(centR_all,1)
nn_Rall_Gs(i) = min(pdist2(centR_all(i,:),centGs));
nn_Rrand_Gs(i) = min(pdist2(centRall_rand(i,:),centGs_rand));
%nn_Pall_Gall(i) = min(pdist2(centP_all(i,:),centG_all));
%nn_Prand_Gall(i) = min(pdist2(centPall_rand(i,:),centG_all));
end
%
disp('startP2')
parfor i=1:size(centRs,1)
   nn_Rs_Gs(i) = min(pdist2(centRs(i,:),centGs));
 %  nn_Ps_Gall(i) = min(pdist2(centPs(i,:),centG_all));
end
disp('startP3')
parfor i=1:size(centRns,1)
   nn_Rns_Gs(i) = min(pdist2(centRns(i,:),centGs));
 %  nn_Pns_Gall(i) = min(pdist2(centPns(i,:),centG_all));
end  
% 
save([base_folder 'nearest_neightbor_pairing_gw10pw10.mat'],'nn_*')
%% Get the plots. 
binsl = 0:0.05:5;
[hy1, hx2] = hist(log10(nn_Gall_Rs),binsl);
[hy2, hx2] = hist(log10(nn_Gs_Rs),binsl);
[hy3, hx2] = hist(log10(nn_Gns_Rs),binsl);
[hy4, hx2] = hist(log10(nn_Grand_Rs),binsl);
%[hy1a, hx2] = hist(log10(nn_Gall_Pall),binsl);
%[hy2a, hx2] = hist(log10(nn_Gs_Pall),binsl);
%[hy3a, hx2] = hist(log10(nn_Gns_Pall),binsl);
%[hy4a, hx2] = hist(log10(nn_Grand_Pall),binsl);

[hy5, hx2] = hist(log10(nn_Rall_Gs),binsl);
[hy6, hx2] = hist(log10(nn_Rs_Gs),binsl);
[hy7, hx2] = hist(log10(nn_Rns_Gs),binsl);
[hy8, hx2] = hist(log10(nn_Rrand_Gs),binsl);
%[hy5a, hx2] = hist(log10(nn_Pall_Gall),binsl);
%[hy6a, hx2] = hist(log10(nn_Ps_Gall),binsl);
%[hy7a, hx2] = hist(log10(nn_Pns_Gall),binsl);
%[hy8a, hx2] = hist(log10(nn_Prand_Gall),binsl);
%
figure;
plot(hx2,hy1,'k'); alpha(0.5); hold on
plot(hx2,hy2,'g'); alpha(0.5); hold on
plot(hx2,hy3,'r'); alpha(0.5); hold on
plot(hx2,hy4*0.5,'b'); alpha(0.5)

savefig([base_folder 'nnGsplit_Rs_counts_with_rand_and_all.fig'])
%%
%figure;
%plot(hx2,hy1a,'k'); alpha(0.5); hold on
%plot(hx2,hy2a,'g'); alpha(0.5); hold on
%plot(hx2,hy3a,'r'); alpha(0.5); hold on
%plot(hx2,hy4a*0.4,'b'); alpha(0.5)
%savefig([base_folder 'nnGsplit_Pall_counts_with_rand_and_all.fig'])

figure;
plot(hx2,hy5,'k'); alpha(0.5); hold on
plot(hx2,hy6,'g'); alpha(0.5); hold on
plot(hx2,hy7,'r'); alpha(0.5); hold on
plot(hx2,hy8*0.5,'b'); alpha(0.5)

savefig([base_folder 'nnRsplit_Gs_counts_with_rand_and_all.fig'])
%
%figure;
%plot(hx2,hy5a,'k'); alpha(0.5); hold on
%plot(hx2,hy6a,'g'); alpha(0.5); hold on
%plot(hx2,hy7a,'r'); alpha(0.5); hold on
%plot(hx2,hy8a*0.4,'b'); alpha(0.5)
%savefig([base_folder 'nnPsplit_Gall_counts_with_rand_and_all.fig'])

%% No longer necessary. 
% ft = fittype( 'gauss2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% %
% clear hyall_use hys_use hyns_use
% hyall_use = cat(1,hy1,hy5);
% hys_use = cat(1,hy2,hy6);
% hyns_use = cat(1,hy3,hy7);
% for i=1:2
%     disp(i)
% [xData, yData] = prepareCurveData( hx2, hyall_use(i,:));
% ft = fittype( 'gauss2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts. Lower = [0 1.8 0 0 2.4 0];
% opts.StartPoint = [12849 2.2 0.3 i^2*15803.883942715 2.8 0.4];
% opts. Upper = [Inf 2.3 0.4 Inf 3.5 0.5];
% 
% [fr, gof] = fit( xData, yData, ft, opts );
% confintall = confint(fr);
% %
% %figure( 'Name', 'untitled fit 1' );
% %subplot(1,3,1)
% %h = plot( fr, xData, yData );
% clear xData yData
% [xData, yData] = prepareCurveData( hx2, hys_use(i,:));
% ft = fittype( 'gauss2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts. Lower = [0 fr.b1 fr.c1 0 fr.b2 fr.c2];
% opts.StartPoint = [hys_use(i,round(fr.b1*10)) fr.b1 fr.c1 hys_use(i,round(fr.b2*10)) fr.b2 fr.c2];
% opts. Upper = [Inf fr.b1*1.001 fr.c1*1.001 Inf fr.b2*1.001 fr.c2*1.001];
% [fitresults, gofs] = fit( xData, yData, ft, opts );
% %subplot(1,3,2)
% %h = plot( fitresults, xData, yData );
% 
% 
% [xData2, yData2] = prepareCurveData( hx2, hyns_use(i,:));
% clear ft opts
% ft = fittype( 'gauss2' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts. Lower = [0 fr.b1 fr.c1 0 fr.b2 fr.c2];
% opts.StartPoint = [hyns_use(i,round(fr.b1*10)) fr.b1 fr.c1 hyns_use(i,round(fr.b2*10)) fr.b2 fr.c2];
% opts. Upper = [Inf fr.b1*1.001 fr.c1*1.001 Inf fr.b2*1.001 fr.c2*1.001];
% [fitresultns, gofns] = fit( xData2, yData2, ft, opts );
% %subplot(1,3,3)
% %h = plot( fitresultns, xData2, yData2 );
% goodfit(i,1) = gof.adjrsquare;
% goodfit(i,2) = gofs.adjrsquare;
% goodfit(i,3) = gofns.adjrsquare;
% %
% sa1(i) = fitresults.a1;
% sb1(i) = fitresults.b1;
% sc1(i) = fitresults.c1;
% sa2(i) = fitresults.a2;
% sb2(i) = fitresults.b2;
% sc2(i) = fitresults.c2;
% nsa1(i) = fitresultns.a1;
% nsb1(i) = fitresultns.b1;
% nsc1(i) = fitresultns.c1;
% nsa2(i) = fitresultns.a2;
% nsb2(i) = fitresultns.b2;
% nsc2(i) = fitresultns.c2;
% syngeph(i) =  fitresults.a1*fitresults.c1 *10*sqrt(2*pi);
% paired_s(i) = fitresults.a1*fitresults.c1 *10*sqrt(2*pi)/...
%     ( fitresults.a1*fitresults.c1 *10*sqrt(2*pi)+ ...
%     fitresults.a2*fitresults.c2 *10*sqrt(2*pi));
% unpaired_s(i) = fitresults.a2*fitresults.c2 *10*sqrt(2*pi)/...
%     ( fitresults.a1*fitresults.c1 *10*sqrt(2*pi)+ ...
%     fitresults.a2*fitresults.c2 *10*sqrt(2*pi));
% paired_ns(i) = fitresultns.a1*fitresultns.c1 *10*sqrt(2*pi)/...
%     ( fitresultns.a1*fitresultns.c1 *10*sqrt(2*pi)+ ...
%     fitresultns.a2*fitresultns.c2 *10*sqrt(2*pi));
% unpaired_ns(i) = fitresultns.a2*fitresultns.c2 *10*sqrt(2*pi)/...
%     ( fitresultns.a1*fitresultns.c1 *10*sqrt(2*pi)+ ...
%     fitresultns.a2*fitresultns.c2 *10*sqrt(2*pi));
% end
% %%
% goodfit
% paired_s
% unpaired_s
% paired_ns
% unpaired_ns
% %%
% figure;
% NumTrials=5; % Number of trail fits per sample
% PeakShape=1; % Gaussian=1, Lorentzian=2, Logistic=3, Pearson=4
% % "extra" determines the shape of the Pearson function. When extra=1
% % the shape is Lorentzian; when extra is large (>20), the shape 
% % approaches Gaussian. At small values (<1), the shape is a pointy "cusp".
% extra=2; 
% center=5; % x-value in center of fitted range
% window =2000;
% NumPeaks = 2;
% %startvector=[100 50 300 100 ];
% startvector=[2 1 2.7 1];
% 
% [FitResultsG,LowestError]=peakfit([hx2' (hy5)'],center,window,NumPeaks,PeakShape,extra,NumTrials,startvector); 

%% Find paired G clusters. 
%Change 'cutoffg' to get better-looking heatmap. 
load([base_folder 'add_to_statsGw10_edges.mat'],'tintsG_p140');
mints_g70s = (([tintsG_p140])./[volumeGs]');
%
val1w = log10(mints_g70s +1)';
val2w = log10(nn_Gs_Rs);  

Xn=70; Yn=80; Xrange=[min(val1w) max(val1w)]; Yrange=[min(val2w) max(val2w)];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ;
%
val1w = log10(mints_g70s +1)';
val2w = log10(nn_Gs_Rs);  
figure; H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange); close;
cutoffg = 60; 
H1 = H; H1(H1>cutoffg)=cutoffg;
figure; 
pcolor(X,Y,H1)

%% OPtics searching. 
k=200;
dataall = (cat(2,val1w,val2w));
datause = dataall(randi(numel(dataall(:,1)),[5000 1]),:);
[RDg,CDg,orderg]=optics(zscore(datause),k);
figure; plot(RDg(orderg))
% figure; plot(CDg(orderg))
%% get dbscan clusters from optics threshold and plot (geph)
%The 'Eps' value should be slightly lower than the first visable peak of
%the figure generated by "Optics searching". 

%It will output a figure with a line split the two group of clusters. 
Eps = 0.22; clustID = 1; classg = zeros(numel(RDg),1);
for i = 1:numel(RDg)
    if RDg(orderg(i))>Eps
        if CDg(orderg(i))<=Eps
            clustID = clustID + 1; classg(orderg(i)) = clustID;
        end
    else
        classg(orderg(i)) = clustID;
    end
end
max(classg)
figure; 
plot(datause(classg==3,1),datause(classg==3,2),'g.'); hold all
plot(datause(classg==2,1),datause(classg==2,2),'r.')
%
figure; pcolor(X,Y,H1); hold on
%ezcontour(@(x,y)pdf(gm,[x y]),[0 180],[3 8]);
%sp(:,1) = (V(1,1)*val1wh + V(2,1)*val2wh)<=-40;
cls = ClassificationDiscriminant.fit(datause,classg);
%Plot the classification boundaries.
K = cls.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
L = cls.Coeffs(2,3).Linear;% boundary between the second and third classes                 
clear f
% Plot the curve K + [x,y]*L  = 0.
xval = [0 200];
yval = -((L(1)/L(2))*xval + K/L(2));
h2 = line(xval, yval);
set(h2,'Color','r','LineWidth',2)
savefig([base_folder 'storm_gs_ps_mints_nn_shell_2d_hist.fig'])

%% Generate output for G. 
%Note that sometimes the line has a negative slop. If that happens the
%following line should be: 
%pairedg_idx = -((L(1)/L(2))*val1w + K/L(2)) < val2w;
pairedg_idx = -((L(1)/L(2))*val1w + K/L(2)) > val2w;
numel(find(pairedg_idx))/numel((pairedg_idx))

%
save([base_folder 'pairing_index_ps_gs_withedges.mat'],'pairedg_idx')

%
load([outpath 'statsG2w10_edges_plus.mat']);
statsGwater_ss = statsGwater(pairedg_idx);
statsGwater_sn = statsGwater(~pairedg_idx);
%
clear tintsG_p140 mints_g140s
load([base_folder 'add_to_statsRw10_edges.mat'],'tintsG_p140');
mints_R70s = (([tintsG_p140])./[volumeRs]');
% Everyhing here is same as above, but for R clsuters. 
Xn=70; Yn=80; Xrange=[-0.2 3]; Yrange=[1.2 4.0];
Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ;
%
val3w = log10(mints_R70s +1)';
val4w = log10(nn_Rs_Gs);  
figure; H = hist2d(cat(2,val3w,val4w),Xn,Yn,Xrange,Yrange); close;
cutoffg = 100; 
H1 = H; H1(H1>cutoffg)=cutoffg;
figure; 
pcolor(X,Y,H1)
%%
k=200;
dataall = (cat(2,val3w,val4w));
datause = dataall(randi(numel(dataall(:,1)),[5000 1]),:);
[RDg,CDg,orderg]=optics(zscore(datause),k);
figure; plot(RDg(orderg))
% figure; plot(CDg(orderg))
%% get dbscan clusters from optics threshold and plot (geph)
Eps = 0.91; clustID = 1; classg = zeros(numel(RDg),1);
for i = 1:numel(RDg)
    if RDg(orderg(i))>Eps
        if CDg(orderg(i))<=Eps
            clustID = clustID + 1; classg(orderg(i)) = clustID;
        end
    else
        classg(orderg(i)) = clustID;
    end
end
max(classg)
figure; 
plot(datause(classg==3,1),datause(classg==3,2),'g.'); hold all
plot(datause(classg==2,1),datause(classg==2,2),'r.')
%
figure; pcolor(X,Y,H1); hold on
%ezcontour(@(x,y)pdf(gm,[x y]),[0 180],[3 8]);
%sp(:,1) = (V(1,1)*val1wh + V(2,1)*val2wh)<=-40;
cls = ClassificationDiscriminant.fit(datause,classg);
%Plot the classification boundaries.
K = cls.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
L = cls.Coeffs(2,3).Linear;% boundary between the second and third classes                 
clear f
% Plot the curve K + [x,y]*L  = 0.
xval = [0 200];
yval = -((L(1)/L(2))*xval + K/L(2));
h2 = line(xval, yval);
set(h2,'Color','r','LineWidth',2)
savefig([base_folder 'storm_rs_gs_mints_nn_shell_2d_hist.fig'])
%%
pairedg_idx = -((L(1)/L(2))*val3w + K/L(2)) > val4w;
disp(numel(find(pairedg_idx))/numel((pairedg_idx)));
save([base_folder 'pairing_index_rs_gs_withedges.mat'],'pairedg_idx')
%
load([outpath 'statsR2w10_edges_plus.mat']);
statsRwater_ss = statsGwater(pairedg_idx);
statsRwater_sn = statsGwater(~pairedg_idx);
%save([base_folder 'R_paired.mat'],'statsRwater_ss','statsRwater_sn','-v7.3');

%%
%From the 0.3% saturated clusters, get corresponded information from 0.1%
%saturated images. 

%If 0.3% saturated images were used from the beginning, then only use the 
% save out part of this block. 
Storm_path = [exp_folder 'analysis\elastic_align_rescale\storm_merged\'];
files = [dir([Storm_path '*.tif'])];
infos = imfinfo([Storm_path files(1,1).name]);
num_images = numel(files);
Storm_images_R = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
Storm_images_G = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');

parfor i = 1:num_images
    Storm_image = imread([Storm_path sprintf('%03d',i) '.tif']);
    Storm_images_R(:,:,i) = Storm_image(:,:,1);
    Storm_images_G(:,:,i) = Storm_image(:,:,2);
end
%
statsRwater_ss_01 = [];
for i = 1:numel(statsRwater_ss)
    stats_temp = stats_copy(statsRwater_ss(i),Storm_images_R);
    statsRwater_ss_01 = cat(2,statsRwater_ss_01,stats_temp);
end
statsGwater_ss_01 = [];
for i = 1:numel(statsGwater_ss)
    stats_temp = stats_copy(statsGwater_ss(i),Storm_images_G);
    statsGwater_ss_01 = cat(2,statsGwater_ss_01,stats_temp);
end
statsRwater_sn_01 = [];
for i = 1:numel(statsRwater_sn)
    stats_temp = stats_copy(statsRwater_sn(i),Storm_images_R);
    statsRwater_sn_01 = cat(2,statsRwater_sn_01,stats_temp);
end
statsGwater_sn_01 = [];
for i = 1:numel(statsGwater_sn)
    stats_temp = stats_copy(statsGwater_sn(i),Storm_images_G);
    statsGwater_sn_01 = cat(2,statsGwater_sn_01,stats_temp);
end

clear statsGwater_ss statsRwater_ss statsGwater_sn statsRwater_sn
statsGwater_ss = statsGwater_ss_01;
statsRwater_ss = statsRwater_ss_01;
statsGwater_sn = statsGwater_sn_01;
statsRwater_sn = statsRwater_sn_01;
%
save([outpath 'G_paired_2.mat'],'statsGwater_ss','statsGwater_sn')
save([outpath 'R_paired_2.mat'],'statsRwater_ss','statsRwater_sn')
disp("data saved")
%% Quantification. 
save_path = 'Y:\Chenghang\ET33_Tigre\20230504_1\analysis\Result\Quantification\';
temp = [statsGwater_ss.Volume1_0];
temp = temp';
save([save_path 'statsGwater_ss.Area.txt'],'temp','-ascii','-double');
temp = [statsGwater_ss.TintsG];
temp = temp';
save([save_path 'statsGwater_ss.TintsG.txt'],'temp','-ascii','-double');
temp = [statsRwater_ss.Volume1_0];
temp = temp';
save([save_path 'statsRwater_ss.Area.txt'],'temp','-ascii','-double');
temp = [statsRwater_ss.TintsG];
temp = temp';
save([save_path 'statsRwater_ss.TintsG.txt'],'temp','-ascii','-double');
temp = [statsGwater_sn.Volume1_0];
temp = temp';
save([save_path 'statsGwater_sn.Area.txt'],'temp','-ascii','-double');
temp = [statsGwater_sn.TintsG];
temp = temp';
save([save_path 'statsGwater_sn.TintsG.txt'],'temp','-ascii','-double');
temp = [statsRwater_sn.Volume1_0];
temp = temp';
save([save_path 'statsRwater_sn.Area.txt'],'temp','-ascii','-double');
temp = [statsRwater_sn.TintsG];
temp = temp';
save([save_path 'statsRwater_sn.TintsG.txt'],'temp','-ascii','-double');
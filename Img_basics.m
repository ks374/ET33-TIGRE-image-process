classdef Img_basics
    %IMG_BASICS 
    
    properties
        Inpath char;
        Img_height double;
        Img_width double;
        num_images int16;
    end
    
    methods
        function obj = Img_basics(Inpath)
            test = split(Inpath,'\');
            if ~isempty(test{end})
                Inpath = [Inpath,'\'];
            end
            test = split(Inpath,'\');
            if string(test{end-1}) ~= "analysis"
                disp('Warning: Make sure the input directory end with .\analysis\');
            end
            obj.Inpath = Inpath;

            temp_path = [Inpath,'elastic_align\for_align\'];
            files = [dir([temp_path '*.png']) dir([temp_path '*.tif'])];
            obj.num_images = numel(files);
            info = imfinfo([temp_path files(1,1).name]);
            obj.Img_height = info.Height;
            obj.Img_width = info.Width;
        end
        
        function img_stack = Img_reader(obj,folder_path,channel_id)
            %channel_id: 0: single channel image; 1-3: single channel from
            %RGB; 4: multichannel RGB. 
            if channel_id ~= 4
                img_stack = zeros(obj.Img_height,obj.Img_width,obj.num_images,'uint8');
            else
                img_stack = zeros(obj.Img_height,obj.Img_width,3,obj.num_images,'uint8');
            end
            if channel_id == 0
                parfor i = 1:obj.num_images
                    Img_temp = imread([folder_path sprintf('%03d',i) '.tif']);
                    img_stack(:,:,i) = Img_temp;
                end
            elseif channel_id < 4
                parfor i = 1:obj.num_images
                    Img_temp = imread([folder_path sprintf('%03d',i) '.tif']);
                    img_stack(:,:,i) = Img_temp(:,:,channel_id);
                end
            else
                parfor i = 1:obj.num_images
                    Img_temp = imread([folder_path sprintf('%03d',i) '.tif']);
                    img_stack(:,:,:,i) = Img_temp;
                end
            end
        end

        function img_stack = load_stats(obj,stats)
            img_stack = zeros(obj.Img_height,obj.Img_width,obj.num_images,'uint8');
            disp('Loading stats to images...')
            for i = 1:numel(stats)
                PixelList = stats(i).PixelList;
                PixelValues = stats(i).PixelValues;
                for j = 1:size(PixelList,1)
                    x = PixelList(j,2);
                    y = PixelList(j,1);
                    z = PixelList(j,3);
                    img_stack(x,y,z) = PixelValues(j);
                end
            end
        end

        function stats_A = add_to(obj,outpath,stats_A,stats_B,save_mat,out_string,radius)
            %Search signal intensity from A to B. 
            %stats_* can be structure or string. 
            %Voxel by default = [15.5,15.5,70], hard-coded. 
            %save_mat: whether or not to save the .mat file.

            if class(stats_A) == "string"
                stats_A = obj.(stats_A);
            end
            if class(stats_B) == "string"
                stats_B = obj.(stats_B);
            end
            Img_H = obj.Img_height;
            Img_W = obj.Img_width;
            Img_N = obj.num_images;
            clear BP BG bg2
            disp('allocating arrays')
            if class(stats_B) == "struct"
                BP2 = zeros(Img_H, Img_W, Img_N,'uint8');
                for i = 1:numel(stats_B)
                    BP2(stats_B(i).PixelIdxList)=stats_B(i).PixelValues;
                end
            else
                BP2 = stats_B;
            end
            namelist = ["tints_p"+string(radius),"volume_p"+string(radius),"area_p"+string(radius) ,...
                "WeightedCentroid_p"+string(radius),"centG_p"+string(radius),...
                "tintsG_p"+string(radius),"volumeG_p"+string(radius)];
            eval(namelist(5)+"=zeros(numel(stats_A),3);");
            eval(namelist(6)+"=zeros(1,numel(stats_A));");
            eval(namelist(7)+"=zeros(1,numel(stats_A));");
            for i=numel(stats_A):-1:1
                stats_A(i).(namelist(1)) = [];
                stats_A(i).(namelist(2)) = [];
                stats_A(i).(namelist(3)) = [];
                stats_A(i).(namelist(4)) = [];
                minpix = min(stats_A(i).PixelList);  maxpix = max(stats_A(i).PixelList);
                min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
                max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
                if min1 < 1; min1=1; end 
                if min2 < 1; min2=1; end
                if min3 < 1; min3=1; end
                if max1 > Img_W; max1=Img_W; end
                if max2 > Img_H; max2=Img_H; end
                if max3 > Img_N; max3=Img_N; end      
                % BPp(i).mat = BP(min2:max2,min1:max1,min3:max3); 
                BP2p(i).mat = BP2(min2:max2,min1:max1,min3:max3); 
            end
            
            disp('computing arrays...');
            parfor jj=1:numel(stats_A)
                %disp(jj);
                voxel = [15.5,15.5,70];
                namelist_par = ["tints_p"+string(radius),"volume_p"+string(radius),"area_p"+string(radius) ,...
                 "WeightedCentroid_p"+string(radius)];
                minpix = min(stats_A(jj).PixelList);  maxpix = max(stats_A(jj).PixelList);
                min1 = minpix(1)-30; min2 = minpix(2)-30; min3 = minpix(3)-6;
                max1 = maxpix(1)+30; max2 = maxpix(2)+30; max3 = maxpix(3)+6;
                if min1 < 1; min1=1; end 
                if min2 < 1; min2=1; end
                if min3 < 1; min3=1; end
                if max1 > Img_W; max1=Img_W; end
                if max2 > Img_H; max2=Img_H; end
                if max3 > Img_N; max3=Img_N; end        
             
                size1 = max1-min1 + 1; size2 = max2-min2 + 1; size3 = max3-min3 + 1;
                curr2 = false(size2, size1, size3);
             
                for j=1: numel(stats_A(jj).PixelList(:,1))
                    curr2(stats_A(jj).PixelList(j,2)-min2+1, ...
                       stats_A(jj).PixelList(j,1)-min1+1, ...
                       stats_A(jj).PixelList(j,3)-min3+1)= 1;  
                end
                % curr1a = BPp(jj).mat;
                curr1b = BP2p(jj).mat;
                Dg = bwdistsc1(curr2,[voxel(1),voxel(2),voxel(3)]);
 
                size2= radius;
                curr3b = Dg<=(size2); curr4b = or(curr2,curr3b);
                [y,x,z] = ind2sub(size(curr4b),find(curr4b));
                PixelList_m140 = [x,y,z];
                RPixelValues_m140 = zeros(numel(PixelList_m140(:,1)),1);
                for kk=1:numel(PixelList_m140(:,1))
                   RPixelValues_m140(kk) = curr1b(PixelList_m140(kk,2),...
                       PixelList_m140(kk,1), PixelList_m140(kk,3));
                end
 
                stats_A(jj).(namelist_par(1)) = sum(RPixelValues_m140);
                stats_A(jj).(namelist_par(3)) = numel(RPixelValues_m140);
                stats_A(jj).(namelist_par(2)) = numel(find(RPixelValues_m140>0));
                stats_A(jj).(namelist_par(4))(1) = sum([PixelList_m140(:,1 )].*...
                        double(RPixelValues_m140))/(sum(RPixelValues_m140) );
                stats_A(jj).(namelist_par(4))(2) = sum([PixelList_m140(:,2 )].*...
                        double(RPixelValues_m140))/(sum(RPixelValues_m140) );
                stats_A(jj).(namelist_par(4))(3) =  sum([PixelList_m140(:, 3)].*...
                        double(RPixelValues_m140))/(sum(RPixelValues_m140));
                
                data_sum_1(jj,:) = stats_A(jj).(namelist_par(4)); %#ok<PFOUS>
                data_sum_2(jj) = stats_A(jj).(namelist_par(1)); %#ok<PFOUS>
                data_sum_3(jj) = stats_A(jj).(namelist_par(2)); %#ok<PFOUS>
            end 
            eval(namelist(5)+"=data_sum_1;");
            eval(namelist(6)+"=data_sum_2;");
            eval(namelist(7)+"=data_sum_3;");
            if save_mat==true
                save([outpath out_string{1} '.mat'],namelist(5),namelist(6),namelist(7),'-v7.3');
                save([outpath out_string{2} '.mat'],'stats_A','-v7.3');
            end
        end

        function pairedg_idx = matcher(obj,outpath,stats_A,stats_B,radius,simple_match_logical)
            % Calculate distance matrxi from the add_to method and idnetify
            % matching. 
            % Note: this method won't calculate the non-selected clusters.
            %simple_match_logical: if true, the match will only be
            %performed based on tints_p* and no distance. 
            if class(stats_A) == "string"
                stats_A = obj.(stats_A);
            end
            
            tints_string = "tints_p"+string(radius);
            tints_A = [stats_A.(tints_string)];
            if simple_match_logical == false
                if class(stats_B) == "string"
                    stats_B = obj.(stats_B);
                end
                disp('Calculating distance matrix');
                cents_A = Get_cents_list(stats_A);
                cents_B = Get_cents_list(stats_B);
                parfor i=1:size(cents_A,1)
                   nn_Gs_Rs(i) = min(pdist2(cents_A(i,:),cents_B));
                end
                volumeAs = [stats_A.Volume1_0];
                mints_g70s = tints_A./volumeAs;
                val1w = log10(mints_g70s +1)';
                val2w = log10(nn_Gs_Rs)';  
    
                Xn=70; Yn=80; Xrange=[min(val1w) max(val1w)]; Yrange=[min(val2w) max(val2w)];
                Xlo = Xrange(1) ; Xhi = Xrange(2) ; Ylo = Yrange(1) ; Yhi = Yrange(2) ; 
                X = linspace(Xlo,Xhi,Xn)' ; Y = linspace(Ylo,Yhi,Yn)' ;
                figure; H = hist2d(cat(2,val1w,val2w),Xn,Yn,Xrange,Yrange); close;
                cutoffg = 10; 
                H1 = H; H1(H1>cutoffg)=cutoffg;
                figure; 
                pcolor(X,Y,H1)
    
                k=100;
                dataall = (cat(2,val1w,val2w));
                datause = dataall(randi(numel(dataall(:,1)),[5000 1]),:);
                [RDg,CDg,orderg]=optics(zscore(datause),k);
                % figure; plot(RDg(orderg))
                
                y_values = RDg(orderg);
                counter = 0;
                for i = 2:numel(y_values)
                    if y_values(i) < (y_values(i-1)-0.05)
                        counter = counter+1;
                    end
                    if counter == 2
                        Eps = y_values(i);
                        break;
                    end
                end
                clustID = 1; classg = zeros(numel(RDg),1);
                for i = 1:numel(RDg)
                    if RDg(orderg(i))>Eps
                        if CDg(orderg(i))<=Eps
                            clustID = clustID + 1; classg(orderg(i)) = clustID;
                        end
                    else
                        classg(orderg(i)) = clustID;
                    end
                end
                if max(classg) < 3
                    disp("Max number of cluster error, check the Optics results");
                end
                figure; 
                plot(datause(classg==3,1),datause(classg==3,2),'g.'); hold all
                plot(datause(classg==2,1),datause(classg==2,2),'r.')
                
                figure; pcolor(X,Y,H1); hold on
                cls = ClassificationDiscriminant.fit(datause,classg);
                K = cls.Coeffs(2,3).Const; % First retrieve the coefficients for the linear
                L = cls.Coeffs(2,3).Linear;% boundary between the second and third classes                 
                clear f
                % Plot the curve K + [x,y]*L  = 0.
                xval = [0 200];
                yval = -((L(1)/L(2))*xval + K/L(2));
                h2 = line(xval, yval);
                set(h2,'Color','r','LineWidth',2)
                savefig([outpath 'storm_gs_ps_mints_nn_shell_2d_hist.fig'])
    
                if -L(1)/L(2) > 0
                    pairedg_idx = -((L(1)/L(2))*val1w + K/L(2)) > val2w;
                else
                    pairedg_idx = -((L(1)/L(2))*val1w + K/L(2)) < val2w;
                end
                disp(numel(find(pairedg_idx))/numel((pairedg_idx)));
            else
                pairedg_idx = tints_A > 0;
            end
            save([outpath 'pairing_index.mat'],'pairedg_idx');
        end

        function H = stats_checker(obj,varargin)
            %stats_checker(obj,block_size,n_seq=1:10,stats_A = "statsGwater_ss",stats_B,stats_C)
            if nargin < 4
                disp("Need more input!");
                return;
            else
                block_size = varargin{1};
                n_seq = varargin{2};
                stats_A = varargin{3};
                stats_A = obj.(stats_A);
                if nargin == 5
                    stats_B = varargin{4};
                    stats_B = obj.(stats_B);
                elseif nargin == 6
                    stats_B = varargin{4};
                    stats_B = obj.(stats_B);
                    stats_C = varargin{5};
                    stats_C = obj.(stats_C);
                end
            end
            
            for i = 1:numel(n_seq)
                cur_stats_A = stats_A(n_seq(i));
                Pixlist = cur_stats_A.PixelList;
                PixValues = cur_stats_A.PixelValues;
                z_start = min(Pixlist(:,3));
                z_end = max(Pixlist(:,3));
                Thickness = z_end - z_start + 1;
                WC = cur_stats_A.WeightedCentroid;
                WC = [WC(2),WC(1),WC(3)];
                box_start = [floor(WC(1)-block_size),floor(WC(2)-block_size)];
                box_end = [floor(WC(1)+block_size),floor(WC(2)+block_size)];
                if box_start(1)<0; box_start(1) = 1; end
                if box_start(2)<0; box_start(2) = 1; end
                if box_end(1)>obj.Img_height; box_end(1)=obj.Img_height; end
                if box_end(2)>obj.Img_width; box_end(2)=obj.Img_width; end
                if nargin == 6
                    Img_block_R = stats_C(box_start(1):box_end(1),box_start(2):box_end(2),z_start:z_end);
                    Img_block_R = max(Img_block_R,[],3);
                else
                    Img_block_R = zeros(box_end(1)-box_start(1)+1,box_end(2)-box_start(2)+1,'uint8');
                end
                if nargin >= 5
                    Img_block_B = stats_B(box_start(1):box_end(1),box_start(2):box_end(2),z_start:z_end);
                    Img_block_B = max(Img_block_B,[],3);
                else
                    Img_block_B = zeros(box_end(1)-box_start(1)+1,box_end(2)-box_start(2)+1,'uint8');
                end
                Img_block_G = zeros(box_end(1)-box_start(1)+1,box_end(2)-box_start(2)+1,Thickness,'uint8');
                for j = 1:numel(PixValues)
                    Img_block_G(Pixlist(j,2)-box_start(1)+1,...
                        Pixlist(j,1)-box_start(2)+1,...
                        Pixlist(j,3)-z_start+1) = PixValues(j);
                end
                Img_block_G = max(Img_block_G,[],3);

                Projection = cat(3,Img_block_R,Img_block_G,Img_block_B);
                H = figure;
                imshow(Projection);
            end
        end

        %WIP: 
        function cluster_restore(obj)
            %Restore the 0.3% saturation for cluster tints property. Also
            %recalculate the volume and weighted centroid. 
        end
        function saver(obj)
            %save the variable in .mat file. 
        end
        function data_writer(obj)
            %write essential properties of the clusters, including: 
            %volume; tints
        end


    end
end

% %%
% %From the 0.3% saturated clusters, get corresponded information from 0.1%
% %saturated images. 
% 
% %If 0.3% saturated images were used from the beginning, then only use the 
% % save out part of this block. 
% Storm_path = [exp_folder 'analysis\elastic_align_rescale\storm_merged\'];
% files = [dir([Storm_path '*.tif'])];
% infos = imfinfo([Storm_path files(1,1).name]);
% num_images = numel(files);
% Storm_images_R = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
% Storm_images_G = zeros(ceil((infos(1,1).Height)),ceil((infos(1,1).Width)),num_images,'uint8');
% 
% parfor i = 1:num_images
%     Storm_image = imread([Storm_path sprintf('%03d',i) '.tif']);
%     Storm_images_R(:,:,i) = Storm_image(:,:,1);
%     Storm_images_G(:,:,i) = Storm_image(:,:,2);
% end
% %
% statsRwater_ss_01 = [];
% for i = 1:numel(statsRwater_ss)
%     stats_temp = stats_copy(statsRwater_ss(i),Storm_images_R);
%     statsRwater_ss_01 = cat(2,statsRwater_ss_01,stats_temp);
% end
% statsGwater_ss_01 = [];
% for i = 1:numel(statsGwater_ss)
%     stats_temp = stats_copy(statsGwater_ss(i),Storm_images_G);
%     statsGwater_ss_01 = cat(2,statsGwater_ss_01,stats_temp);
% end
% statsRwater_sn_01 = [];
% for i = 1:numel(statsRwater_sn)
%     stats_temp = stats_copy(statsRwater_sn(i),Storm_images_R);
%     statsRwater_sn_01 = cat(2,statsRwater_sn_01,stats_temp);
% end
% statsGwater_sn_01 = [];
% for i = 1:numel(statsGwater_sn)
%     stats_temp = stats_copy(statsGwater_sn(i),Storm_images_G);
%     statsGwater_sn_01 = cat(2,statsGwater_sn_01,stats_temp);
% end
% 
% clear statsGwater_ss statsRwater_ss statsGwater_sn statsRwater_sn
% statsGwater_ss = statsGwater_ss_01;
% statsRwater_ss = statsRwater_ss_01;
% statsGwater_sn = statsGwater_sn_01;
% statsRwater_sn = statsRwater_sn_01;
% %
% save([outpath 'G_paired_2.mat'],'statsGwater_ss','statsGwater_sn')
% save([outpath 'R_paired_2.mat'],'statsRwater_ss','statsRwater_sn')
% disp("data saved")
% %% Quantification. 
% save_path = 'Y:\Chenghang\ET33_Tigre\20230504_1\analysis\Result\Quantification\';
% temp = [statsGwater_ss.Volume1_0];
% temp = temp';
% save([save_path 'statsGwater_ss.Area.txt'],'temp','-ascii','-double');
% temp = [statsGwater_ss.TintsG];
% temp = temp';
% save([save_path 'statsGwater_ss.TintsG.txt'],'temp','-ascii','-double');
% temp = [statsRwater_ss.Volume1_0];
% temp = temp';
% save([save_path 'statsRwater_ss.Area.txt'],'temp','-ascii','-double');
% temp = [statsRwater_ss.TintsG];
% temp = temp';
% save([save_path 'statsRwater_ss.TintsG.txt'],'temp','-ascii','-double');
% temp = [statsGwater_sn.Volume1_0];
% temp = temp';
% save([save_path 'statsGwater_sn.Area.txt'],'temp','-ascii','-double');
% temp = [statsGwater_sn.TintsG];
% temp = temp';
% save([save_path 'statsGwater_sn.TintsG.txt'],'temp','-ascii','-double');
% temp = [statsRwater_sn.Volume1_0];
% temp = temp';
% save([save_path 'statsRwater_sn.Area.txt'],'temp','-ascii','-double');
% temp = [statsRwater_sn.TintsG];
% temp = temp';
% save([save_path 'statsRwater_sn.TintsG.txt'],'temp','-ascii','-double');

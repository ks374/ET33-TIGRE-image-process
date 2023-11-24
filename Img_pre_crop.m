classdef Img_pre_crop < Img_basics
    %Crop the image before further processing. 
    %   
    
    properties
        Folderlist cell;

        ROI double;
        target_image uint8;
    end
    
    methods
        function obj = Img_pre_crop(Inpath,folderlist)
            obj = obj@Img_basics(Inpath);
            if nargin < 2
                obj.Folderlist = {'conv_488','for_align','conv_merged','storm_merged'};
            else
                obj.Folderlist = folderlist;
            end

        end
        
        function obj= Img_reader(obj,folder_id,channel_id)
            %Read images from folder specified by the folder_id and
            %channel_id into the target image stack
            %Channel_id: 0 single color image; 1-3 RGB image but read one
            %color;4 RGC image read all images. 
            folder = [obj.Inpath,'elastic_align\',obj.Folderlist{folder_id},'\'];
            obj.target_image = [];
            obj.target_image = Img_reader@Img_basics(obj,folder,channel_id);
            
        end

        function obj= Img_reader_ds(obj,folder_id,channel_id)
            %Read images from folder specified by the folder_id and
            %channel_id(varargin) into the target image stack
            folder = [obj.Inpath,'elastic_align\',obj.Folderlist{folder_id},'_ds\'];
            img_stack = zeros(ceil(obj.Img_height/10),ceil(obj.Img_width/10),obj.num_images,'uint8');
            for i = 1:obj.num_images
                Img_temp = imread([folder sprintf('%03d',i) '.tif']);
                if channel_id == 0
                    img_stack(:,:,i) = Img_temp;
                else
                    img_stack(:,:,i) = Img_temp(:,:,channel_id);
                end
            end
            obj.target_image = img_stack;
        end

        function img_proj = find_project(obj)
            img_proj = max(obj.target_image,[],3);
        end

        function ROI = selection(~,img_proj)
            %Draw a rectangle on the image projection and store it in the
            %ROI variable
            % ROI = [Starting_y,Starting_x,length_y,length_x]. 
            [~,ROI] = imcrop(img_proj);
            close;
        end

        function obj = ROI_overlap(obj,varargin)
            %Find the overlapping region of selected ROIs
            ROI = zeros(numel(varargin),4);
            n_input = numel(varargin);
            for i = 1:n_input
                ROI(i,:) = varargin{i};
                ROI(i,3) = ROI(i,1) + ROI(i,3);
                ROI(i,4) = ROI(i,2) + ROI(i,4);
            end
            ROI(1,1) = max(ROI(:,1));
            ROI(2,1) = max(ROI(:,2));
            ROI(3,1) = min(ROI(:,3));
            ROI(4,1) = min(ROI(:,4));
            ROI = ROI(1,:);
            ROI(1,3) = ROI(1,3) - ROI(1,1);
            ROI(1,4) = ROI(1,4) - ROI(1,2);

            obj.ROI = ROI*10;
        end

        function batch_save(obj)
            %Read the images, crop based on ROI, and save it to the target
            %folder.
            % Note: ROIs are downscaled ROI. 
            outpath_base = [obj.Inpath 'elastic_align_rescale\'];
            if ~exist(outpath_base,'dir')
                mkdir(outpath_base);
            else
                rmdir(outpath_base,'s');
                mkdir(outpath_base);
            end
            ROI = obj.ROI;
            save([outpath_base,'Cropper.mat'],'ROI');
            for i = 1:numel(obj.Folderlist)
                outpath = [outpath_base obj.Folderlist{i} '\'];
                outpath_ds = [outpath_base obj.Folderlist{i} '_ds\'];
                mkdir(outpath);
                mkdir(outpath_ds);
                splitted_foldername = split(obj.Folderlist{i},'_');
                if string(splitted_foldername{end}) == "merged"
                    channel_id = 4;
                    obj = obj.Img_reader(i,channel_id);
                    for j = 1:obj.num_images
                        out_image = imcrop(obj.target_image(:,:,:,j),obj.ROI);
                        imwrite(out_image,[outpath sprintf('%03d',j) '.tif']);
                        imwrite(imresize(out_image,0.1),[outpath_ds sprintf('%03d',j) '.tif']);
                    end
                else
                    channel_id = 0;
                    obj = obj.Img_reader(i,channel_id);
                    for j = 1:obj.num_images
                        out_image = imcrop(obj.target_image(:,:,j),obj.ROI);
                        imwrite(out_image,[outpath sprintf('%03d',j) '.tif']);
                        imwrite(imresize(out_image,0.1),[outpath_ds sprintf('%03d',j) '.tif']);
                    end
                end
            end
        end
    end
end


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
            for i = 1:numel(stats)
                disp('Loading stats to images...')
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

        function add_to(obj,A,B)
            %An abstract add_to function
            pass;
        end

    end
end


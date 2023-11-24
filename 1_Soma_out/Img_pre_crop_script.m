exp_folder = 'Y:\Chenghang\ET33_Tigre\20230504_1\analysis\';
Folderlist = {'conv_488','for_align','conv_merged','storm_merged'};

Cropper = Img_pre_crop(exp_folder,Folderlist);
Cropper = Cropper.Img_reader_ds(3,1);%Red for soma
Img_proj = Cropper.find_project;
ROI_1 = Cropper.selection(Img_proj);
Cropper = Cropper.Img_reader_ds(4,3);%Blue for Axons
Img_proj = Cropper.find_project;
ROI_2 = Cropper.selection(Img_proj);

Cropper = Cropper.ROI_overlap(ROI_1,ROI_2);

Cropper.batch_save;
clear;clc
Input_path = 'Y:\Chenghang\ET33_Tigre\20230504_1\analysis\';
outpath = [Input_path 'Result\4_S_on_axon\'];
%%
%Manually find VGLuT2 clusters that are centeredin in convetional axonal
%staining: 
block_size = 100;
n_seq = 1:100;
stats_A = "statsGwater_ss";
stats_B = "Img_axon";
stats_C = "Img_bassoon";
AS_conv_0 = Axon_searcher(Input_path,true);
AS_conv_0.stats_checker(block_size,n_seq,stats_A,stats_B,stats_C);
n_seq = 6;
AS_STORM = Axon_searcher(Input_path,false);
AS_conv_0.stats_checker_both(block_size,n_seq,stats_A,stats_B,stats_C,AS_STORM);
%%
outfile = [outpath 'Pair_Result.xlsx'];
search_radius = 0;
AS_conv_0 = Axon_searcher(Input_path,true);
AS_conv_0 = AS_conv_0.add_to(outpath,"statsGwater_ss","Img_axon",0,{},search_radius);
AS_conv_0 = AS_conv_0.matcher(outpath,"statsGwater_ss","Img_axon",search_radius,true);
AS_conv_0.match_writter(outfile,"statsGwater_ss",search_radius,"Conv");
search_radius = 70;
AS_conv_70 = AS_conv_0.add_to(outpath,"statsGwater_ss","Img_axon",0,{},search_radius);
AS_conv_70 = AS_conv_70.matcher(outpath,"statsGwater_ss","Img_axon",search_radius,true);
AS_conv_70.match_writter(outfile,"statsGwater_ss",search_radius,"Conv");
search_radius = 140;
AS_conv_140 = AS_conv_0.add_to(outpath,"statsGwater_ss","Img_axon",0,{},search_radius);
AS_conv_140 = AS_conv_140.matcher(outpath,"statsGwater_ss","Img_axon",search_radius,true);
AS_conv_140.match_writter(outfile,"statsGwater_ss",search_radius,"Conv");
search_radius = 210;
AS_conv_210 = AS_conv_0.add_to(outpath,"statsGwater_ss","Img_axon",0,{},search_radius);
AS_conv_210 = AS_conv_210.matcher(outpath,"statsGwater_ss","Img_axon",search_radius,true);
AS_conv_210.match_writter(outfile,"statsGwater_ss",search_radius,"Conv");
search_radius = 280;
AS_conv_280 = AS_conv_0.add_to(outpath,"statsGwater_ss","Img_axon",0,{},search_radius);
AS_conv_280 = AS_conv_280.matcher(outpath,"statsGwater_ss","Img_axon",search_radius,true);
AS_conv_280.match_writter(outfile,"statsGwater_ss",search_radius,"Conv");
search_radius = 350;
AS_conv_350 = AS_conv_0.add_to(outpath,"statsGwater_ss","Img_axon",0,{},search_radius);
AS_conv_350 = AS_conv_350.matcher(outpath,"statsGwater_ss","Img_axon",search_radius,true);
AS_conv_350.match_writter(outfile,"statsGwater_ss",search_radius,"Conv");
%%
%VGLuT2 clusters that are matched with Bassoon is more likely to be on
%conventional axon staining? 
ss_idx = [AS_conv.statsGwater_ss.paired_idx];
sn_idx = [AS_conv.statsGwater_sn.paired_idx];
disp(sum(ss_idx)/numel(ss_idx));
disp(sum(sn_idx)/numel(sn_idx));
%%
AS_Storm_0 = Axon_searcher(Input_path,false);
%
outfile = [outpath 'Pair_Storm_0.xlsx'];
search_radius = 0;
AS_Storm_0 = AS_Storm_0.add_to(outpath,"statsGwater_ss","statsVa2s",0,{},search_radius);
AS_Storm_0 = AS_Storm_0.matcher(outpath,"statsGwater_ss","statsVa2s",search_radius,false);
AS_Storm_0.match_writter(outfile,"statsGwater_ss",search_radius);
outfile = [outpath 'Pair_Storm_70.xlsx'];
search_radius = 70;
AS_Storm_70 = AS_Storm_0.add_to(outpath,"statsGwater_ss","statsVa2s",0,{},search_radius);
AS_Storm_70 = AS_Storm_70.matcher(outpath,"statsGwater_ss","statsVa2s",search_radius,false);
AS_Storm_70.match_writter(outfile,"statsGwater_ss",search_radius);
outfile = [outpath 'Pair_Storm_140.xlsx'];
search_radius = 140;
AS_Storm_140 = AS_Storm_0.add_to(outpath,"statsGwater_ss","statsVa2s",0,{},search_radius);
AS_Storm_140 = AS_Storm_140.matcher(outpath,"statsGwater_ss","statsVa2s",search_radius,false);
AS_Storm_140.match_writter(outfile,"statsGwater_ss",search_radius);
outfile = [outpath 'Pair_Storm_210.xlsx'];
search_radius = 210;
AS_Storm_210 = AS_Storm_0.add_to(outpath,"statsGwater_ss","statsVa2s",0,{},search_radius);
AS_Storm_210 = AS_Storm_210.matcher(outpath,"statsGwater_ss","statsVa2s",search_radius,false);
AS_Storm_210.match_writter(outfile,"statsGwater_ss",search_radius);
outfile = [outpath 'Pair_Storm_280.xlsx'];
search_radius = 280;
AS_Storm_280 = AS_Storm_0.add_to(outpath,"statsGwater_ss","statsVa2s",0,{},search_radius);
AS_Storm_280 = AS_Storm_280.matcher(outpath,"statsGwater_ss","statsVa2s",search_radius,false);
AS_Storm_280.match_writter(outfile,"statsGwater_ss",search_radius);
outfile = [outpath 'Pair_Storm_350.xlsx'];
search_radius = 350;
AS_Storm_350 = AS_Storm_0.add_to(outpath,"statsGwater_ss","statsVa2s",0,{},search_radius);
AS_Storm_350 = AS_Storm_350.matcher(outpath,"statsGwater_ss","statsVa2s",search_radius,false);
AS_Storm_350.match_writter(outfile,"statsGwater_ss",search_radius);
%%
AS_conv.stats_checker_both(100,1:10,"statsGwater_ss","Img_axon","Img_bassoon",AS_Storm);

%WIP: write for non-simple conveitonal image matching

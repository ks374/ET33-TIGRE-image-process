classdef Axon_searcher < Img_basics
    %AXON_SEARCHER Identify whether a synaptic cluster is on Storm or
    %conventional axonal staining. 
    
    properties
        Inpath char;
        
        statsGwater_ss struct;
        statsGwater_sn struct;
        statsVa2s struct;
        Img_axon_storm uint8;
        Img_axon_conv uint8;
    end
    
    methods
        function obj = Axon_searcher(Inpath)
            obj = obj@Img_basics(Inpath);
            
            Result_path = [Inpath 'Result\'];
            load([Result_path 'G_paired_2.mat'],"statsGwater_ss","statsGwater_sn");
            obj.statsGwater_ss = statsGwater_ss;
            obj.statsGwater_sn = statsGwater_sn;
            load([Inpath 'Result\3_Axon\statsV2sw10.mat'],"statsGa2s");
            obj.statsVa2s = statsGa2s;

            Conv_axon_path = [Inpath 'elastic_align_rescale\conv_488\'];
            obj.Img_axon_conv = obj.Img_reader(Conv_axon_path,0);
            obj.Img_axon_storm = obj.load_stats(obj,statsVa2s);
        end
        
        function outputArg = method1(obj,inputArg)
            
        end
    end
end


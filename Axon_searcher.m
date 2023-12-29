classdef Axon_searcher < Img_basics
    %AXON_SEARCHER Identify whether a synaptic cluster is on Storm or
    %conventional axonal staining. 
    
    properties
        statsGwater_ss struct;
        statsGwater_sn struct;
        statsVa2s struct;
        Img_axon uint8;
        Img_bassoon uint8;
    end
    
    methods
        function obj = Axon_searcher(Inpath,logical_conv)
            obj = obj@Img_basics(Inpath);
            
            Result_path = [Inpath 'Result\'];
            load([Result_path 'G_paired_2.mat'],"statsGwater_ss","statsGwater_sn");
            
            obj.statsGwater_ss = statsGwater_ss;
            obj.statsGwater_sn = statsGwater_sn;
            load([Inpath 'Result\3_Axon\statsV2sw10.mat'],"statsGa2s");
            obj.statsVa2s = statsGa2s;

            Conv_axon_path = [Inpath 'elastic_align\conv_488_thred\'];
            if logical_conv == true
                obj.Img_axon = obj.Img_reader(Conv_axon_path,0);
            else
                obj.Img_axon = obj.load_stats(obj.statsVa2s);
            end

            load([Result_path 'R_paired_2.mat'],"statsRwater_ss");
            obj.Img_bassoon = obj.load_stats(statsRwater_ss);
        end
        
        function obj = add_to(obj,outpath,stats_A,stats_B,save_mat,out_string,radius)
            stats_temp = add_to@Img_basics(obj,outpath,stats_A,stats_B,save_mat,out_string,radius);
            obj.(stats_A) = stats_temp;
        end
        
        function obj = matcher(obj,outpath,stats_A,stats_B,radius,simple_match_logical)
            paired_idx = matcher@Img_basics(obj,outpath,stats_A,stats_B,radius,simple_match_logical);
            for i = 1:numel(obj.(stats_A))
                obj.(stats_A)(i).paired_idx = paired_idx(i);
            end
        end

        function stats_checker(obj,block_size,n_seq,stats_A,stats_B,stats_C)
            %stats_checker(obj,block_size,n_seq=1:10,stats_A = "statsGwater_ss",stats_B,stats_C)
            stats_temp = obj.(stats_A);
            for i = 1:numel(n_seq)
                disp("Checking stats# "+string(n_seq(i))+":");
                %disp("Stats paired_idx = " + string(stats_temp(n_seq(i)).paired_idx));
                H = stats_checker@Img_basics(obj,block_size,n_seq(i),stats_A,stats_B,stats_C);
                %saveas(H,)
                %close H;
            end
        end

        function stats_checker_both(obj,block_size,n_seq,stats_A,stats_B,stats_C,AS_another)
            %check conventional and STORM at the same time
            for i = 1:numel(n_seq)
                disp("Checking stats# "+string(n_seq(i))+":");
                %stats_temp = obj.(stats_A);
                %disp("Stats Conv paired_idx = " + string(stats_temp(n_seq(i)).paired_idx));
                obj.stats_checker(block_size,n_seq(i),stats_A,stats_B,stats_C);
                %stats_temp = AS_another.(stats_A);
                %disp("Stats Storm paired_idx = " + string(stats_temp(n_seq(i)).paired_idx));
                AS_another.stats_checker(block_size,n_seq(i),stats_A,stats_B,stats_C);
                disp("Press any key to continue...");pause;close all;
            end

        end

        function match_writter(obj,outfile,stats_A,search_radius)
            %Write the pairing info to a csv file. 
            stats_temp = obj.(stats_A);
            Headline = ["Cluster_ID","search_radius","Pair_id"];
            writematrix(Headline,outfile);
            for i = 1:numel(stats_temp)
                Textline = [string(i),string(search_radius),string(stats_temp(i).paired_idx)];
                writematrix(Textline,outfile,'WriteMode','append');
            end
            disp("Done file writting");
        end

    end
end


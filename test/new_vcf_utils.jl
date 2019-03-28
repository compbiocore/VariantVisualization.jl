#=
using RoguePkg

julia> using VariantVisualization
INFO: Loading HttpServer methods...
WARNING: Method definition ==(Base.Nullable{S}, Base.Nullable{T}) in module Base at nullable.jl:238 overwritten in module NullableArrays at /Users/George/.julia/test/v0.6/NullableArrays/src/operators.jl:99.

julia> Pkg.test(pkg_for"VariantVisualization")

=#

@testset "VCFUtils" begin

vcf_filename = "test_files/test_4X_191.vcf"
#vcf_filename_with_chr = "test_files/test_with_chr.vcf"

reader = VCF.Reader(open(vcf_filename, "r"))
#reader_with_chr = VCF.Reader(open(vcf_filename_with_chr, "r"))
sample_names = get_sample_names(reader)

dp_limit = 50

@testset "clean_column1!" begin
    df = Matrix(["X" 1 2; "Y" 2 3; 2 4 1])
    clean_column1!(df)
    @test df[1,1] == "23"
    @test df[2,1] == "24"
    @test df[3,1] == 2
end

#functions for variant filters

@testset "io_genomic_range_vcf_filter" begin
sub = io_genomic_range_vcf_filter("chr4:0-400000000",vcf_filename)
@test typeof(sub) == Array{Any,1}
@test size(sub,1) == 1012
#println("io_genomic_range_vcf_filter type is $(typeof(sub))")
#println("io_genomic_range_vcf_filter size is $(size(sub,1))")
end


@testset "filters_with_siglist" begin

    @testset "load_siglist" begin
    sig_list=load_siglist("test_files/positions_list.csv")

    @testset "pass_genomic_range_siglist_filter" begin
    sig_list=load_siglist("test_files/positions_list.csv")
    sub = pass_genomic_range_siglist_filter(vcf_filename,sig_list,"chr4:0-5000000000")

    @test (typeof(sub[1])) == GeneticVariation.VCF.Record
    @test (length(sub)) ==  5

    end

            @testset "io_sig_list_vcf_filter" begin

            sub = VariantVisualization.io_sig_list_vcf_filter(sig_list,vcf_filename)
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  11
            end

            @testset "pass_siglist_filter" begin
            sub = pass_siglist_filter(vcf_filename, sig_list)
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  10
            end

            @testset "genomic_range_siglist_filter" begin
            sub = genomic_range_siglist_filter(vcf_filename,sig_list,"chr4:0-400000000")
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  5
            end

    end
end


@testset "io_pass_filter" begin
    sub = io_pass_filter(vcf_filename)
    @test (typeof(sub[1])) == GeneticVariation.VCF.Record
    @test (length(sub)) ==  1164
end

@testset "pass_genomic_range_filter" begin
    reader = VCF.Reader(open(vcf_filename, "r"))
    sub = pass_genomic_range_filter(reader,"chr4:0-400000000",vcf_filename)
    @test (typeof(sub[1])) == GeneticVariation.VCF.Record
    @test (length(sub)) ==  856
end

#functions for converting vcf record array to numerical array
@testset "combined_all_genotype_array_functions" begin

sub = io_pass_filter(vcf_filename)
number_rows = size(sub,1)

gt_num_array,gt_chromosome_labels=combined_all_genotype_array_functions(sub)
#println("combined_all_genotype_array_functions gt array is type: $(typeof(gt_num_array))")
#println("combined_all_genotype_array_functions gt array is length: $(length(gt_num_array))")
#println("combined_all_genotype_array_functions gt_chromosome_labels is typeof: $(typeof(gt_chromosome_labels))")
#println("combined_all_genotype_array_functions gt_chromosome_labels is length: $(length(gt_chromosome_labels))")
@test typeof(gt_num_array) == Array{Int64,2}
@test length(gt_num_array) == 222324
@test typeof(gt_chromosome_labels) == Array{Any,2}
@test length(gt_chromosome_labels) == 2328

    @testset "chromosome_label_generator" begin
    chrom_label_info = VariantVisualization.chromosome_label_generator(gt_chromosome_labels[:,1])
    #println("chromosome_label_generator chrom_label_info is type $(typeof(chrom_label_info))")
    #println("chromosome_label_generator chrom_label_info is length $(length(chrom_label_info))")
    @test typeof(chrom_label_info) == Tuple{Array{String,1},Array{Int64,1},String}
    @test size(chrom_label_info,1) == 3
    end

    @testset "generate_chromosome_positions_for_hover_labels" begin
    chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(gt_chromosome_labels)
    #println("generate_chromosome_positions_for_hover_labels chr_pos_tuple_list is type $(typeof(chr_pos_tuple_list))")
    #println("generate_chromosome_positions_for_hover_labels chr_pos_tuple_list is length $(length(chr_pos_tuple_list))")
    @test typeof(chr_pos_tuple_list) == Array{String,1}
    @test size(chr_pos_tuple_list,1) == 1164
    end

    @testset "generate_genotype_array" begin
    sub = io_pass_filter(vcf_filename)
    genotype_array=generate_genotype_array(sub,"GT")

    #println("generate_genotype_array is $(typeof(genotype_array))")
    #println("generate_genotype_array is $(size(genotype_array,1))")
    @test typeof(genotype_array) == Array{String,2}
    @test size(genotype_array,1) == 1164

    @testset "define_geno_dict" begin
    geno_dict = define_geno_dict()
    #println("define_geno_dict is type is $(typeof(geno_dict))")
    #println("define_geno_dict is length is $(length(geno_dict))")
    @test typeof(geno_dict) == Dict{Any,Any}
    @test length(geno_dict) == 100

    @testset "translate_genotype_to_num_array" begin
    gt_num_array,gt_chromosome_labels=translate_genotype_to_num_array(genotype_array,geno_dict)
    #println("translate_genotype_to_num_array gt_num_array type is $(typeof(gt_num_array))")
    #println("translate_genotype_to_num_array gt_num_array length is $(length(gt_num_array))")
    #println("translate_genotype_to_num_array gt_chromosome_labels typeof is $(typeof(gt_chromosome_labels))")
    #println("translate_genotype_to_num_array gt_chromosome_labels length is $(length(gt_chromosome_labels))")
    @test typeof(gt_num_array) == Array{Int64,2}
    @test length(gt_num_array) == 222324
    @test typeof(gt_chromosome_labels) == Array{String,2}
    @test length(gt_chromosome_labels) == 2328
    end
    end
    end

end

@testset "combined_all_read_depth_array_functions" begin

sub = io_pass_filter(vcf_filename)
dp_num_array,dp_chromosome_labels=combined_all_read_depth_array_functions(sub)

#println("combined_all_read_depth_array_functions dp_num_array type is $(typeof(dp_num_array))")
#println("combined_all_read_depth_array_functions dp_num_array length is $(length(dp_num_array))")
#println("combined_all_read_depth_array_functions gt_chromosome_labels typeof is $(typeof(dp_chromosome_labels))")
#println("combined_all_read_depth_array_functions gt_chromosome_labels length is $(length(dp_chromosome_labels))")
@test typeof(dp_num_array) == Array{Int64,2}
@test length(dp_num_array) == 222324
@test typeof(dp_chromosome_labels) == Array{Any,2}
@test length(dp_chromosome_labels) == 2328

    @testset "get_sample_names" begin
    reader = VCF.Reader(open(vcf_filename, "r"))
    sample_names=get_sample_names(reader)
    #println("get_sample_names")
    #println("get_sample_names sample_names type is $(typeof(sample_names))")
    #println("get_sample_names sample_names length is $(length(sample_names))")
    @test typeof(sample_names) == Array{Symbol,2}
    @test length(sample_names) == 191

            @testset "read_depth_threshhold" begin
            dp_num_array=read_depth_threshhold(dp_num_array)
            #println("read_depth_threshhold dp_num_array is type $(typeof(dp_num_array))")
            #println("read_depth_threshhold dp_num_array is size $(size(dp_num_array,1))")
            @test typeof(dp_num_array) == Array{Int64,2}
            @test size(dp_num_array,1) == 1164
            end

        @testset "avg_dp_samples" begin
        avg_sample_list=avg_dp_samples(dp_num_array)
        #println("avg_dp_samples avg_sample_list type is $(typeof(avg_sample_list))")
        #println("avg_dp_samples avg_sample_list length is $(length(avg_sample_list))")
        @test typeof(avg_sample_list) == Array{Float64,1}
        @test size(avg_sample_list,1) == 191

        @testset "list_sample_names_low_dp" begin
        list=list_sample_names_low_dp(avg_sample_list,sample_names)
        #println("list_sample_names_low_dp list type is $(typeof(list))")
        #println("list_sample_names_low_dp list length is $(length(list))")
        @test typeof(list) == Array{String,1}
        @test size(list,1) == 6
        end

    end

    @testset "avg_dp_variant" begin
    avg_variant_list=avg_dp_variant(dp_num_array)
    #println("avg_dp_variant avg_variant_list type is $(typeof(avg_variant_list))")
    #println("avg_dp_variant avg_variant_list length is $(length(avg_variant_list))")
    @test typeof(avg_variant_list) == Array{Float64,1}
    @test size(avg_variant_list,1) == 1164

        @testset  "list_variant_positions_low_dp" begin
        list=list_variant_positions_low_dp(avg_variant_list,dp_chromosome_labels)
        #println("list_variant_positions_low_dp list type is $(typeof(list))")
        #println("list_variant_positions_low_dp list length is $(length(list))")
        @test typeof(list) == Array{Tuple{Any,Int64},1}
        @test size(list,1) == 33

        end

    end

    @testset "sortcols_by_phenotype_matrix" begin
    vcf,group_label_pack=sortcols_by_phenotype_matrix("test_files/sample_metadata_matrix.csv","case,control", dp_num_array, sample_names)
    #println("sortcols_by_phenotype_matrix vcf type is $(typeof(vcf))")
    #println("sortcols_by_phenotype_matrix vcf size is $(size(vcf,1))")
    #println("sortcols_by_phenotype_matrix group_label_pack type is $(typeof(group_label_pack))")
    #println("sortcols_by_phenotype_matrix group_label_pack size is $(size(group_label_pack,1))")
    @test typeof(vcf) == Array{Int64,2}
    @test size(vcf,1) == 1164
    @test typeof(group_label_pack) == Array{Any,1}
    @test size(group_label_pack,1) == 5

#=
        @testset "find_group_label_indices" begin
        pheno = DelimitedFiles.readdlm("test_files/fixed_pheno_matrix_test.csv", ',')
        row_to_sort_by = find(x -> x == "control,case", pheno)
        row_to_sort_by = row_to_sort_by[1]
        group_label_pack=find_group_label_indices(pheno,"control,case",row_to_sort_by)
        #println("find_group_label_indices group_label_pack type is $(typeof(group_label_pack))")
        #println("find_group_label_indices group_label_pack length is $(length(group_label_pack))")
        @test typeof(group_label_pack) == Array{Any,1}
        @test size(group_label_pack,1) == 5
        end

        =#

        @testset "select_columns" begin
        dp_num_array=select_columns("test_files/select_samples_list.txt", dp_num_array, sample_names)
        #println("select_columns dp_num_array type is $(typeof(dp_num_array))")
        #println("select_columns dp_num_array size is $(size(dp_num_array,1))")
        @test typeof(dp_num_array) == Tuple{Array{Int64,2},Array{Any,1}}
        @test size(dp_num_array,1) == 2
        end

    end

    end
    end

end

#=
using RoguePkg

julia> using ViVa
INFO: Loading HttpServer methods...
WARNING: Method definition ==(Base.Nullable{S}, Base.Nullable{T}) in module Base at nullable.jl:238 overwritten in module NullableArrays at /Users/George/.julia/test/v0.6/NullableArrays/src/operators.jl:99.

julia> Pkg.test(pkg_for"ViVa")

=#

@testset "VCFUtils" begin

vcf_filename = "test_files/test_4X_191.vcf"
vcf_filename_with_chr = "test_files/test_with_chr.vcf"

reader = VCF.Reader(open(vcf_filename, "r"))
reader_with_chr = VCF.Reader(open(vcf_filename_with_chr, "r"))
sample_names = get_sample_names(reader)

@testset "clean_column1!" begin
    df = Matrix(["X" 1 2; "Y" 2 3; 2 4 1])
    clean_column1!(df)
    @test df[1,1] == "23"
    @test df[2,1] == "24"
    @test df[3,1] == 2
end

#functions for variant filters

@testset "io_chromosome_range_vcf_filter" begin
sub = io_chromosome_range_vcf_filter("chr4:0-400000000",reader)
println(sub[1:2])
println(size(sub,2))
end

#=
@testset "filters_with_siglist" begin

    @testset "load_siglist" begin
    sig_list=load_siglist("test_files/significantList_for_proteinstructures.csv")
    println(sig_list[2:1])
    println(size(sig_list,1))

            @testset "clean_column1_siglist!" begin
            clean_column1_siglist!(sig_list)
            println(sig_list[1,2])
            println(size(sig_list,1))
            end

            @testset "io_sig_list_vcf_filter" begin
            sub=io_sig_list_vcf_filter(sig_list,vcf_filename)
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  13
            end

            @testset "pass_chrrange_siglist_filter" begin
            sub = pass_chrrange_siglist_filter(vcf_filename,sig_list,"chr4:0-400000000")
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  12
            end

            @testset "pass_siglist_filter" begin
            sub = pass_siglist_filter(vcf_filename, sig_list)
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  12
            end

            @testset "chrrange_siglist_filter" begin
            sub = chrrange_siglist_filter(vcf_filename,sig_list,"chr4:0-400000000")
            @test (typeof(sub[1])) == GeneticVariation.VCF.Record
            @test (length(sub)) ==  13
            end

    end

    end
=#

@testset "io_pass_filter" begin
    reader = VCF.Reader(open(vcf_filename, "r"))
    sub = io_pass_filter(reader)
    @test (typeof(sub[1])) == GeneticVariation.VCF.Record
    @test (length(sub)) ==  1164
end

@testset "pass_chrrange_filter" begin
    reader = VCF.Reader(open(vcf_filename, "r"))
    sub = pass_chrrange_filter(reader,"chr4:0-400000000")
    @test (typeof(sub[1])) == GeneticVariation.VCF.Record
    @test (length(sub)) ==  856
end

#functions for converting vcf record array to numerical array
@testset "combined_all_genotype_array_functions" begin

reader = VCF.Reader(open(vcf_filename, "r"))
sub = io_pass_filter(reader)

gt_num_array,gt_chromosome_labels=combined_all_genotype_array_functions(sub)
println(typeof(gt_num_array))
println(length(gt_num_array))
println(typeof(gt_chromosome_labels))
println(length(gt_chromosome_labels))

    @testset "generate_genotype_array" begin
    reader = VCF.Reader(open(vcf_filename, "r"))
    sub = io_pass_filter(reader)
    genotype_array=generate_genotype_array(sub,"GT")
    println(typeof(genotype_array))
    println(length(genotype_array))
    println(genotype_array[3:5])

    @testset "define_geno_dict" begin
    geno_dict = define_geno_dict()
    println(typeof(geno_dict))
    println(length(geno_dict))

    @testset "translate_genotype_to_num_array" begin
    gt_num_array,gt_chromosome_labels=translate_genotype_to_num_array(genotype_array,geno_dict)
    println(typeof(gt_num_array))
    println(length(gt_num_array))
    println(typeof(gt_chromosome_labels))
    println(length(gt_chromosome_labels))
    end
    end
    end

end

@testset "combined_all_read_depth_array_functions" begin #inside functions same used in combined_all_genotype_array_functions

reader = VCF.Reader(open(vcf_filename, "r"))
sub = io_pass_filter(reader)
dp_num_array,dp_chromosome_labels=combined_all_read_depth_array_functions(sub)
println(typeof(dp_num_array))
println(length(dp_num_array))
println(typeof(dp_chromosome_labels))
println(length(dp_chromosome_labels))

@testset "get_sample_names" begin
reader = VCF.Reader(open(vcf_filename, "r"))
sample_names=get_sample_names(reader)
println("get_sample_names")
println(typeof(sample_names))
println(length(sample_names))

@testset "avg_dp_samples" begin
avg_sample_list=avg_dp_samples(dp_num_array)
println("avg_sample_list is $avg_sample_list")

@testset "list_sample_names_low_dp" begin
list=list_sample_names_low_dp(avg_sample_list,sample_names)
println(list)
end

end



@testset "avg_dp_variant" begin
avg_variant_list=avg_dp_variant(dp_num_array)
println("avg_dp_variant is $avg_variant_list")
end

@testset "sortcols_by_phenotype_matrix" begin
vcf,group_label_pack=sortcols_by_phenotype_matrix("test_files/sample_phenotype_matrix.csv","control,case", dp_num_array, sample_names)
println(typeof(vcf))
println(size(vcf,1))
println(typeof(group_label_pack))
println(length(group_label_pack))

    @testset "find_group_label_indices" begin
    pheno = readdlm("test_files/sample_phenotype_matrix.csv", ',')
    row_to_sort_by = find(x -> x == "control,case", pheno)
    row_to_sort_by = row_to_sort_by[1]
    group_label_pack=find_group_label_indices(pheno,"control,case",row_to_sort_by)
    println(typeof(group_label_pack))
    println(length(group_label_pack))
    end

    @testset "select_columns" begin
    dp_num_array=select_columns("test_files/select_samples_list.txt", dp_num_array, sample_names)
    println(typeof(dp_num_array))
    println(length(dp_num_array))
    end

end
end
end

"""
    list_variant_positions_low_dp(variant_avg_list::Array{Float64,2},chrom_labels)
finds variant positions that have an average read depth of under 15 across all patients
"""
function list_variant_positions_low_dp(variant_avg_list::Array{Float64,1},chrom_labels)

    low_dp_index_list = Array{Int64}(0)

        for item = 1:length(variant_avg_list)
            if variant_avg_list[item] < 15
                push!(low_dp_index_list,item)
            end
        end

    low_dp_positions = Array{Tuple{Int64,Int64}}(0)

        for i in low_dp_index_list
            chrom_position = chrom_labels[i,1],chrom_labels[i,2]
            push!(low_dp_positions,chrom_position)
        end

        return low_dp_positions
end

#functions for producing objects for plot functions

"""
    read_depth_threshhold(dp_array::Array{Int64,2})
sets ceiling for read depth values at dp = 100. All dp over 100 are set to 100 to visualize read depth values between 0 < dp > 100 in better definition
"""
function read_depth_threshhold(dp_array::Array{Int64,2})

    dp_array[dp_array[:,:].>100].=100

    return dp_array
end

"""
    save_numerical_array(num_array::Matrix{Any},sample_names,chr_labels)
save numerical array with chr labels and sample ids to working directory
"""
function save_numerical_array(num_array,sample_names,chr_labels)

      #samplenames=sample_names
      #samplenames=Matrix(samplenames)

      headings = hcat("chr","position")
      sample_names = hcat(headings,sample_names)

      chr_labeled_array_for_plotly=hcat(chr_labels, num_array)
      labeled_value_matrix_withsamplenames= vcat(sample_names,chr_labeled_array_for_plotly)

      writedlm("AC_gatk406_eh_PASS_withheader_value_matrix_.txt", labeled_value_matrix_withsamplenames, "\t")
end

"""
    chromosome_label_generator(chromosome_labels::Array{String,2})
Returns vector of chr labels and indices to mark chromosomes in plotly heatmap
Specifically, saves indexes and chrom labels in vectors to pass into heatmap function to ticvals and tictext respectively
Input is either gt_chromosome_labels or dp_chromosome_labels from translate_gt/dp_to_num_array()
"""

function chromosome_label_generator(chromosome_labels::Array{Any,1})
    chrom_label_indices = findfirst.(map(a -> (y -> isequal(a, y)), unique(chromosome_labels)), [chromosome_labels])
    chrom_labels = unique(chromosome_labels)
    chrom_labels = [string(i) for i in chrom_labels]

    if length(chrom_labels) > 1
        for item=2:(length(chrom_labels))

            ratio=((chrom_label_indices[item])-(chrom_label_indices[item-1]))/(length(chromosome_labels))
            println(ratio)

            if ratio < 0.2
                font_size = "8"
                println("font size is $font_size")
                return chrom_labels,chrom_label_indices,font_size
            else
                font_size = "10"
                println("font size is $font_size")

                return chrom_labels,chrom_label_indices,font_size
            end
        end
    else

        font_size = "10"
        return chrom_labels,chrom_label_indices,font_size
    end
end

"""
    checkfor_outputdirectory(path::String)
Checks to see if output directory exists already. If it doesn't, it creates the new directory to write output files to.
"""
function checkfor_outputdirectory(path::String)
    if isdir(path) == true
    else
            mkdir(path)
    end
end

"""
    generate_chromosome_positions_for_hover_labels(chr_labels::Array{Any,2})
creates tuple of genomic locations to set as tick labels. This is automatically store chromosome positions in hover labels. However tick labels are set to hidden with showticklabels=false so they will not crowd the y axis.
"""
function generate_chromosome_positions_for_hover_labels(chr_labels::Array{Any,2})

returnXY_column1!(chr_labels) #not working yet
#println(chr_labels)

chr_pos_tuple_list=Array{Tuple}(0)

    for row = 1:size(chr_labels,1)

        chr=chr_labels[row,1]
        pos=chr_labels[row,2]
        chr_pos_tuple=chr,pos
        push!(chr_pos_tuple_list,chr_pos_tuple)
    end

    return chr_pos_tuple_list
end

"""
    returnXY_column1!(chr_label_vector)
Replace String "23","24","25" with "X","Y","M" in chromosome label vector used for plot labels
"""
@testset "returnXY_column1!" begin

end

"""
    sort_genotype_array(genotype_array)
sorts genotype array for GT or DP by chromosomal location
"""
@testset sort_genotype_array(genotype_array) begin

end


=#
end

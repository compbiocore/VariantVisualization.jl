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

@testset "io_chromosome_range_vcf_filter" begin
sub = io_chromosome_range_vcf_filter("chr4:0-400000000",reader)
println(sub[1:2])
println(size(sub,2))
end

@testset "io_sig_list_vcf_filter" begin

    @testset "load_siglist" begin
        sig_list=load_siglist("test_files/significantList_for_proteinstructures.csv")
        println(sig_list[2:1])
        println(size(sig_list,1))

            @testset "clean_column1_siglist!" begin
            clean_column1_siglist!(sig_list)
            println(sig_list[1,2])
            println(size(sig_list,1))
            end

        sub=io_sig_list_vcf_filter(sig_list,vcf_filename)
        println(sub[1,5])

        @testset "pass_chrrange_siglist_filter" begin
        sub = pass_chrrange_siglist_filter(vcf_filename,sig_list,"chr4:0-400000000")
        println(sub[1,5])

        @testset pass_siglist_filter begin
        sub = pass_siglist_filter(vcf_filename,sig_list)
        end

        @testset "chrrange_siglist_filter" begin
        sub = chrrange_siglist_filter(vcf_filename,sig_list,"chr4:0-400000000")
        end


        end

    end
end

@testset "io_pass_filter" begin
    sub = io_pass_filter(reader)
    println(sub[2,1])
end

@testset "pass_chrrange_filter" begin
    sub = pass_chrrange_filter(reader,"chr4:0-400000000")
end






#=
#functions for variant filters


#functions for converting vcf record array to numerical array

"""
    create_chr_dict()
creates dict for use in combined_all_genotype_array_functions() for removing 'chr' from chromosome labels to allow sorting variant records by chromosome position.
"""
@testset create_chr_dict() begin

end

"""
    combined_all_genotype_array_functions(sub)
convert sub from variant filters to gt_num_array and gt_chromosome_labels for plot functions.
"""
@testset combined_all_genotype_array_functions(sub) begin

end

"""
    combined_all_read_depth_array_functions(sub)
convert sub from variant filters to dp_num_array and dp_chromosome_labels for plot functions.
"""
@testset combined_all_read_depth_array_functions(sub) begin

end

"""
    generate_genotype_array(record_sub::Array{Any,1},genotype_field::String)
Returns numerical array of genotype values (either genotype or read_depth values) which are translated by another function into num_array
Where genotype_field is either GT or DP to visualize genotype or read_depth
"""
@testset generate_genotype_array(record_sub::Array{Any,1},y) begin

end

"""
    define_geno_dict()
returns dictionary of values for use in replace_genotype_with_vals()
"""
@testset define_geno_dict() begin

end

"""
    translate_genotype_to_num_array(genotype_array,geno_dict)
returns a tuple of num_array for plotting, and chromosome labels for plotting as label bar.
Translates array of genotype values to numerical array of categorical values.
Genotype values are converted to categorical values. No_call=0, 0/0=400, heterozygous_variant=600, homozygous_variant=800
"""
@testset translate_genotype_to_num_array(genotype_array,geno_dict) begin

end

"""
    translate_readdepth_strings_to_num_array(read_depth_array::Array{Any,2})
Returns array of read_depth as int for plotting and average calculation.
By default, read depth values over 100 are replaced with 100 to improve heatmap visualization (see read_depth_threshhold() ).
Where read_depth_array is output of generate_genotype_array() for DP option
returns a tuple of num_array type Int for average calculation and plotting, and chromosome labels for plotting as label bar
"""
@testset translate_readdepth_strings_to_num_array(read_depth_array::Array{Any,2}) begin

end


#functions for sample filters

"""
    get_sample_names(reader)
returns sample ids of vcf file as a vector of symbols for naming columns of num_array dataframe object for column filter functions
"""
@testset get_sample_names(reader) begin

end

"""
    find_group_label_indices(pheno)
find indices and determines names for group 1 and group 2 labels on plots. finds index of center of each sample group to place tick mark and label.
"""
@testset find_group_label_indices(pheno,trait_to_group_by,row_to_sort_by) begin

end

"""
    sortcols_by_phenotype_matrix(pheno_matrix_filename::String,trait_to_group_by::String,num_array::Array{Int64,2}, sample_names::Array{Symbol,2})
group samples by a common trait using a user generated key matrix ("phenotype matrix")
"""
@testset sortcols_by_phenotype_matrix(pheno_matrix_filename::String,trait_to_group_by::String, num_array::Array{Int64,2}, sample_names::Array{Symbol,2}) begin

end

"""
    select_columns(filename_sample_list::AbstractString, num_array::Array{Int64,2}, sample_names::Array{Symbol,2})
returns num_array with columns matching user generated list of sample ids to select for analysis. num_array now has sample ids in first row.
"""
@testset select_columns(filename_sample_list::AbstractString, num_array::Array{Int64,2}, sample_names::Array{Symbol,2}) begin

end


#functions for mathematic analysis
"""
    avg_dp_samples(dp_num_array::Array{Int64,2})
create sample_avg_list vector that lists averages of read depth for each sample for input into avg_sample_dp_line_chart(sample_avg_list)
dp_num_array must contain dp values as Int64 and be without chromosome position columns
"""
@testset avg_dp_samples(dp_num_array::Array{Int64,2}) begin

end


"""
    avg_dp_variant(dp_num_array::Array{Int64,2})
create variant_avg_list vector that lists averages of read depth for each variant for input into avg_variant_dp_line_chart(variant_avg_list)
"""
@testset avg_dp_variant(dp_num_array::Array{Int64,2}) begin

end

"""
    list_sample_names_low_dp(sample_avg_list::Array{Float64,2},sample_names)
returns list of sample ids that have an average read depth of under 15 across all variant positions
"""
@testset list_sample_names_low_dp(sample_avg_list::Array{Float64,1},sample_names) begin

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

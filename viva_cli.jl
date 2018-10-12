#=
read_depth heatmap up to 15 white, up to 30 light blue,
don't output files in same folder
=#

println("loading packages:")
println("ViVa")
using ViVa

println("GeneticVariation")
using GeneticVariation

println("ArgParse")
using ArgParse

println("VCFTools")
using VCFTools


#=
tic()
println("loading packages:")
println("ViVa")
using ViVa
toc()
tic()
println("GeneticVariation")
using GeneticVariation
toc()
tic()
println("ArgParse")
using ArgParse
toc()
tic()
println("VCFTools")
using VCFTools
toc()
=#

tic()

function test_parse_main(ARGS::Vector{String})

    # initialize the settings (the description is for the help screen)
    s = ArgParseSettings(

    description = "ViVa VCF Visualization Tool is a tool for creating publication quality plots of data contained within VCF files. For a complete description of features with examples read the docs here https://github.com/compbiocore/ViVa.jl",
    suppress_warnings = true,
    epilog = "Thank you for using ViVa. Please submit any bugs to https://github.com/compbiocore/ViVa.jl/issues "

    )

    @add_arg_table s begin
        "--vcf_file", "-f"             # vcf filename
        help = "vcf filename in format: file.vcf"
        arg_type = String
        required = true

        "--show_stats"
        help = "show number of records, number samples, etc. in vcf file"
        action = :store_true

        "--output_directory", "-o"
        help =" function checks if directory exists and saves there, if not creates and saves here"
        arg_type = String
        default = "output"

        "--save_format", "-s"               # format to save graphics in
        help = "file format you wish to save graphics as (eg. pdf, html, png)"
        arg_type = String
        required = true

        "--chromosome_range", "-r"
        help = "select rows within a given chromosome range. Provide chromosome range after this flag."
        arg_type = String

        "--pass_filter", "-p"
        help = "select rows within a given chromosome range. Provide chromosome range after this flag in format chr4:20000000-30000000"
        arg_type = String
        action = :store_true

        "--positions_list", "-l"
        help = "select variants matching list of chromosomal positions. Provide filename of text file formatted with two columns in csv format: 1,2000345."
        arg_type = String

        "--group_samples", "-g"
        help = "group samples by common trait using user generated matrix key of traits and sample names following format guidelines in documentation. Provide file name of .csv file"
        nargs = 2
        arg_type = String

        "--select_samples", "-x"
        help = "select samples to include in visualization by providing tab delimited list of sample names (eg. samplenames.txt)"
        arg_type = String

        "--heatmap", "-m"
        help = "genotype field to visualize (eg. genotype, read_depth, or 'genotype,read_depth' to visualize each separately)"
        arg_type = String
        default = "genotype"

        "heatmap_title"
        help = "positional argument. Specify filename for heatmap"
        arg_type = String

        "--avg_dp"
        help = "visualize average read depths as line chart. Options: average sample read depth, average variant read depth, or both. eg. =sample, =variant, =sample,variant"
        arg_type = String

    end

    parsed_args = parse_args(s)# can turn off printing parsed args after development
    #println("Parsed args:")

#activate block to show all argument keys
#=
    for (key,val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end
=#

    return parsed_args

end

parsed_args = test_parse_main(ARGS)

#begin main program

#filter vcf and load matrix

vcf_filename = (parsed_args["vcf_file"])
println("Reading $vcf_filename")
println()

reader = VCF.Reader(open(vcf_filename, "r"))
sample_names = get_sample_names(reader)

ViVa.checkfor_outputdirectory(parsed_args["output_directory"])

if parsed_args["show_stats"] == true

    number_records = nrecords((parsed_args["vcf_file"]))
    number_samples = nsamples((parsed_args["vcf_file"]))

    println("_______________________________________________")
    println()
    println("Summary Statistics of $(parsed_args["vcf_file"])")
    println()
    println("number of records: $number_records")
    println("number of samples: $number_samples")
    println("_______________________________________________")
    println()

end

tic()

if parsed_args["pass_filter"] == true && parsed_args["chromosome_range"] == nothing && parsed_args["positions_list"] == nothing
    sub = ViVa.io_pass_filter(reader)
    number_rows = size(sub,1)
    println("selected $number_rows variants with Filter status: PASS")
    heatmap_input = "pass_filtered"

end

if parsed_args["chromosome_range"] != nothing && parsed_args["pass_filter"] == false && parsed_args["positions_list"] == nothing
    sub = ViVa.io_chromosome_range_vcf_filter(parsed_args["chromosome_range"],reader)
    number_rows = size(sub,1)
    println("selected $number_rows variants within chromosome range: $(parsed_args["chromosome_range"])")
    heatmap_input = "range_filtered"

end

if parsed_args["positions_list"] != nothing && parsed_args["chromosome_range"] == nothing && parsed_args["pass_filter"] == false
    sig_list =  load_siglist(parsed_args["positions_list"])
    sub = ViVa.io_sig_list_vcf_filter(sig_list,reader)
    number_rows = size(sub,1)
    println("selected $number_rows variants that match list of chromosome positions of interest")
    heatmap_input = "positions_filtered"
end

if parsed_args["pass_filter"] == true && parsed_args["chromosome_range"] != nothing && parsed_args["positions_list"] != nothing
    sig_list =  load_siglist(parsed_args["positions_list"])
    sub = ViVa.pass_chrrange_siglist_filter(reader, sig_list, parsed_args["chromosome_range"])
    number_rows = size(sub,1)
    println("selected $number_rows variants with Filter status: PASS, that match list of chromosome positions of interest, and are within chromosome range: $(parsed_args["chromosome_range"])")
end

if parsed_args["pass_filter"] == true && parsed_args["chromosome_range"] != nothing && parsed_args["positions_list"] == nothing
    sub = ViVa.pass_chrrange_filter(reader, parsed_args["chromosome_range"])
    number_rows = size(sub,1)
    println("selected $number_rows variants with Filter status: PASS and are within chromosome range: $(parsed_args["chromosome_range"])")
end

if parsed_args["pass_filter"] == true && parsed_args["chromosome_range"] == nothing && parsed_args["positions_list"] != nothing
    sig_list =  load_siglist(parsed_args["positions_list"])
    sub = ViVa.pass_siglist_filter(reader, sig_list)
    number_rows = size(sub,1)
    println("selected $number_rows variants with Filter status: PASS and that match list of chromosome positions of interest")
end

if parsed_args["pass_filter"] == false && parsed_args["chromosome_range"] == nothing && parsed_args["positions_list"] == nothing

    println("no filters applied. Large vcf files will take a long time to process and heatmap visualizations will lose resolution at this scale unless viewed in interactive html for zooming.")
    println()
    println("Loading VCF file into memory for visualization")

    sub = Array{Any}(0)
    for record in reader
        push!(sub,record)
    end

end

#sub = clean_column1_chr(sub)

if parsed_args["heatmap"] == "genotype"
    gt_num_array,gt_chromosome_labels = combined_all_genotype_array_functions(sub)

    if parsed_args["heatmap_title"] != nothing
        title = parsed_args["heatmap_title"]
    else
        title = "Genotype_$(parsed_args["vcf_file"])"
    end

    chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(gt_chromosome_labels)

    chrom_label_info = ViVa.chromosome_label_generator(gt_chromosome_labels[:,1])

    if length(parsed_args["group_samples"]) == 2

        if parsed_args["select_samples"] != nothing
            println("selecting samples listed in $(parsed_args["select_samples"])")
            gt_num_array = select_columns(parsed_args["select_samples"], gt_num_array, sample_names)
        end

        group_trait_matrix_filename=((parsed_args["group_samples"])[1])
        trait_to_group_by = ((parsed_args["group_samples"])[2])
        println("grouping samples by $trait_to_group_by")
        gt_num_array,group_label_pack = sortcols_by_phenotype_matrix(group_trait_matrix_filename, trait_to_group_by, gt_num_array, sample_names)

        graphic = ViVa.genotype_heatmap_with_groups(gt_num_array,title,chrom_label_info,group_label_pack,sample_names)

    else

        if parsed_args["select_samples"] != nothing
            gt_num_array = select_columns(parsed_args["select_samples"], gt_num_array, sample_names)
        end

        graphic = ViVa.genotype_heatmap2(gt_num_array,title,chrom_label_info,sample_names)
    end

    PlotlyJS.savefig(graphic, joinpath("$(parsed_args["output_directory"])" ,"$title.$(parsed_args["save_format"])"))

end

if parsed_args["heatmap"] == "read_depth"

    dp_num_array,dp_chromosome_labels = combined_all_read_depth_array_functions(sub)

    if parsed_args["heatmap_title"] != nothing
        title = parsed_args["heatmap_title"]
    else
        title = "Read_Depth_$(parsed_args["vcf_file"])"
    end

    chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(dp_chromosome_labels)

    chrom_label_info = ViVa.chromosome_label_generator(dp_chromosome_labels[:,1])

    if length(parsed_args["group_samples"]) == 2

        if parsed_args["select_samples"] != nothing
            println("selecting samples listed in $(parsed_args["select_samples"])")
            dp_num_array = select_columns(parsed_args["select_samples"], dp_num_array, sample_names)
        end

        group_trait_matrix_filename=((parsed_args["group_samples"])[1])
        trait_to_group_by = ((parsed_args["group_samples"])[2])
        println("grouping samples by $trait_to_group_by")

        dp_num_array,group_label_pack = sortcols_by_phenotype_matrix(group_trait_matrix_filename, trait_to_group_by, dp_num_array, sample_names)
        dp_num_array_limited=read_depth_threshhold(dp_num_array)

        graphic = ViVa.dp_heatmap2_with_groups(dp_num_array_limited,title,chrom_label_info,group_label_pack)

    else

        if parsed_args["select_samples"] != nothing
            dp_num_array = select_columns(parsed_args["select_samples"], dp_num_array, sample_names)
        end

        dp_num_array_limited=read_depth_threshhold(dp_num_array)

        graphic = ViVa.dp_heatmap2(dp_num_array_limited, title, chrom_label_info, sample_names)
    end

    PlotlyJS.savefig(graphic, joinpath("$(parsed_args["output_directory"])" ,"$title.$(parsed_args["save_format"])"))

end

#println("_______________________________________________")
println()
println("Finished Filtering. Total time to filter:")
toc()
println("_______________________________________________")
#println()

if parsed_args["avg_dp"] == "sample"

    if isdefined(:dp_num_array) == false
        read_depth_array = ViVa.generate_genotype_array(sub,"DP")
        clean_column1!(read_depth_array)
        read_depth_array=ViVa.sort_genotype_array(read_depth_array)
        dp_num_array,dp_chromosome_labels = translate_readdepth_strings_to_num_array(read_depth_array)
    end

    chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(dp_chromosome_labels)

    avg_list = ViVa.avg_dp_samples(dp_num_array)
    list = ViVa.list_sample_names_low_dp(avg_list, sample_names)
    writedlm(joinpath("$(parsed_args["output_directory"])","Samples_with_low_dp.txt"),list, ",")
    #println("The following samples have read depth of under 15: $list")
    graphic = avg_sample_dp_scatter(avg_list,sample_names)
    PlotlyJS.savefig(graphic, joinpath("$(parsed_args["output_directory"])" ,"Average Sample Read Depth.$(parsed_args["save_format"])"))

elseif parsed_args["avg_dp"] == "variant"
    if isdefined(:dp_num_array) == false
        read_depth_array = ViVa.generate_genotype_array(sub,"DP")
        clean_column1!(read_depth_array)
        read_depth_array=ViVa.sort_genotype_array(read_depth_array)
        dp_num_array,dp_chromosome_labels = translate_readdepth_strings_to_num_array(read_depth_array)
    end

    chr_pos_tuple_list = generate_chromosome_positions_for_hover_labels(dp_chromosome_labels)

    avg_list = ViVa.avg_dp_variant(dp_num_array)
    list = ViVa.list_variant_positions_low_dp(avg_list, dp_chromosome_labels)
    writedlm(joinpath("$(parsed_args["output_directory"])","Variant_positions_with_low_dp.txt"),list,",")
    #println("The following samples have read depth of less than 15: $list")
    graphic = avg_variant_dp_line_chart(avg_list,chr_pos_tuple_list)
    #PlotlyJS.savefig(graphic, "Average Variant Read Depth.$(parsed_args["save_format"])") #make unique save format - default to pdf but on my computer html
    PlotlyJS.savefig(graphic, joinpath("$(parsed_args["output_directory"])" ,"Average Variant Read Depth.$(parsed_args["save_format"])"))
end

close(reader)


#println("_______________________________________________")
println()
println("Finished plotting. Total time to run:")
toc()
#println("_______________________________________________")
println()

#=
println("_______________________________________________")
println()
println("Finished plotting. Total time to run: $(toc()) seconds")
println("_______________________________________________")
println()

=#



#in the morning - write function to convert number array matrix to dataframe for input into column filter functions, get all sample names from io for use here
#add positional argument to save list of positions and samples with dp under 15 as file instead of printint
#add positional arguments to save any num array with labels

#convert records matrix to genotype or read_depth matrix

#if no filter then load vcf and convert to number matrix

#=
outline of command line interface

1) Read / Load vcf
if no variant filters, check size of file and if over threshhold, apply filters (see manual: Variant Filtering)
print warning that vcf file may be too large to store in local memory and needs to be filtered
We don't believe there is a reason to visualize genotype or read depth data at this scale.

2) Variant Filters


3) Convert array of records to genotype field numerical array



4) Sample Selection and Ordering


5) Plotting
    A) array for plotly

    B) create heatmaps for field selections



if


=#

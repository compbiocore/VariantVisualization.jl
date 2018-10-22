#functions for loading vcf, vcf stats, clean vcf, load siglist

"""
    clean_column1!(matrix_with_chr_column)
Replace String "X","Y","M" from chromosome column with 23,24,25 respectively so variants can be sorted by descending chr position for plotting
"""
function clean_column1!(matrix_with_chr_column)

    for i = 1:(size(matrix_with_chr_column, 1))
        matrix_with_chr_column[i, 1] = matrix_with_chr_column[i, 1] == "X" ? "23" : matrix_with_chr_column[i, 1]
        matrix_with_chr_column[i, 1] = matrix_with_chr_column[i, 1] == "Y" ? "24" : matrix_with_chr_column[i, 1]
        matrix_with_chr_column[i, 1] = matrix_with_chr_column[i, 1] == "M" ? "25" : matrix_with_chr_column[i, 1]
    end

end

"""
    clean_column1_chr(matrix_with_chr_column)
Replace String "chr" from chromosome column with "" so chr position is loaded as int and variants can be sorted by descending chr position for plotting
"""
function clean_column1_chr(matrix_with_chr_column)

    for i = 1:(size(matrix_with_chr_column, 1))
        replace(matrix_with_chr_column[i, 1],"chr","")
    end

end

"""
    clean_column1_siglist!(siglist)
    replaces "X","Y","M" with 23,24,25 {Int}
use in load_siglist() because X and Y need to be replaced with Int
"""
function clean_column1_siglist!(siglist)

    for i = 1:(size(siglist, 1))
        siglist[i, 1] = siglist[i, 1] == "X" ? 23 : siglist[i, 1]
        siglist[i, 1] = siglist[i, 1] == "Y" ? 24 : siglist[i, 1]
        siglist[i, 1] = siglist[i, 1] == "M" ? 25 : siglist[i, 1]
    end
end

"""
    returnXY_column1!(chr_label_vector)
Replace String "23","24","25" with "X","Y","M" in chromosome label vector used for plot labels
"""
function returnXY_column1!(chr_label_vector)

    for i = 1:(size(chr_label_vector, 1))
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == "23" ? "X" : chr_label_vector[i, 1]
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == "24" ? "Y" : chr_label_vector[i, 1]
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == "25" ? "M" : chr_label_vector[i, 1]
    end

end

"""
    sort_genotype_array(genotype_array)
sorts genotype array for GT or DP by chromosomal location
"""
function sort_genotype_array(genotype_array)

    data=genotype_array[:,3:size(genotype_array,2)]
    chrom_positions = [parse(Int, i) for i in genotype_array[:,1:2]]
    genotype_array = hcat(chrom_positions,data)

    genotype_array = sortrows(genotype_array, by=x->(x[1],x[2]))

return genotype_array
end

"""
    load_siglist(filename::AbstractString)
where x = filename of significant SNP variant location list in comma delimited format (saved as .csv)
"""
function load_siglist(filename::AbstractString)

siglist_unsorted = readdlm(filename, ',', skipstart=1)
ViVa.clean_column1_siglist!(siglist_unsorted)
siglist = sortrows(siglist_unsorted, by=x->(x[1],x[2]))

if typeof(siglist) == Array{Float64,2}
    siglist = trunc.(Int,siglist)
    return siglist
end

return siglist
end



#functions for variant filters
"""
io_chromosome_range_vcf_filter(chr_range::String, reader::GeneticVariation.VCF.Reader)
create subarray of vcf variant records matching user specified chromosome range in format: (e.g. chr1:0-30000000)
"""
function io_chromosome_range_vcf_filter(chr_range::String,reader::GeneticVariation.VCF.Reader)
       a=split(chr_range,":")
       chrwhole=a[1]
       chrnumber=split(chrwhole,"r")
       string_chr=chrnumber[2]
       chr=String(string_chr)
       range=a[2]
       splitrange=split(range, "-")
       lower_limit=splitrange[1]
       chr_range_low=parse(lower_limit)
       upper_limit=splitrange[2]
       chr_range_high=parse(upper_limit)

       vcf_subarray = Array{Any}(0)

       for record in reader

              if (VCF.chrom(record) == chr) && (chr_range_high > VCF.pos(record) > chr_range_low)
                     push!(vcf_subarray,record)
              end
       end

       return vcf_subarray
end

"""
    io_sig_list_vcf_filter(sig_list,vcf_filename)
returns subarray of variant records matching a list of variant positions returned from load_siglist()
"""
function io_sig_list_vcf_filter(sig_list,vcf_filename)

       vcf_subarray = Array{Any}(0)

       for row= 1:size(sig_list,1)
              dimension = size(sig_list,1)

              chr=(sig_list[row,1])
              pos=(sig_list[row,2])

              reader = VCF.Reader(open(vcf_filename, "r"))

              for record in reader

                     if typeof(VCF.chrom(record)) == String
                            chr = string(chr)

                            if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                                push!(vcf_subarray,record)
                            end

                    else

                            if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                                push!(vcf_subarray,record)

                            end
                    end
              end
       end

       return vcf_subarray
end

#=
function io_sig_list_vcf_filter(sig_list,reader::GeneticVariation.VCF.Reader)

       vcf_subarray = Array{Any}(0)
       println(sig_list)

       for row= 1:size(sig_list,1)
              dimension = size(sig_list,1)

              chr=(sig_list[row,1])
              pos=(sig_list[row,2])

              println("siglist iteration shows chromosome $chr")
              println("siglist iteration shows position $pos")


              for record in reader

#reader = VCF.Reader(open("test_with_chr.vcf", "r"))

                 # println(VCF.chrom(record))
                  println("record iteration shows chromosome $chr")
                  #println(VCF.pos(record))
                  println("record iteration shows pos $pos")

                  println()

                    if typeof(VCF.chrom(record)) == String
                        if typeof(chr) == String
                            chr = String(chr)

                            if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                                println("match!")
                                push!(vcf_subarray,record)
                            end

                        elseif typeof(chr) == Int64
                            chr = string(chr)

                            if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                                    println("match!")
                                    push!(vcf_subarray,record)
                            end
                        end

                    else

                            if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                                println(VCF.chrom(record))
                                println(VCF.pos(record))
                                push!(vcf_subarray,record)

                            end
                    end
              end
       end

       return vcf_subarray
end
=#
"""
    io_pass_filter(reader::GeneticVariation.VCF.Reader)
returns subarray of vcf records including only records with FILTER status = PASS
"""
function io_pass_filter(reader::GeneticVariation.VCF.Reader)

vcf_subarray = Array{Any}(0)

for record in reader

       if VCF.hasfilter(record) && VCF.filter(record) == String["PASS"]
              push!(vcf_subarray,record)
       end
end

return vcf_subarray
end

"""
    pass_chrrange_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)
returns subarray of vcf records with io_pass_filter, io_sig_list_vcf_filter, and io_chromosome_range_vcf_filter applied.
"""
function pass_chrrange_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)

    a=split(chr_range,":")
    chrwhole=a[1]
    chrnumber=split(chrwhole,"r")
    string_chr=chrnumber[2]
    chr=String(string_chr)
    range=a[2]
    splitrange=split(range, "-")
    lower_limit=splitrange[1]
    chr_range_low=parse(lower_limit)
    upper_limit=splitrange[2]
    chr_range_high=parse(upper_limit)

    vcf_subarray = Array{Any}(0)

    for row= 1:size(sig_list,1)
           dimension = size(sig_list,1)

           chr=(sig_list[row,1])
           pos=(sig_list[row,2])

           reader = VCF.Reader(open(vcf_filename, "r"))

           for record in reader

                  if typeof(VCF.chrom(record)) == String
                         chr = string(chr)

                         if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos) && (VCF.hasfilter(record)) && (VCF.filter(record) == String["PASS"]) && ((VCF.chrom(record) == chr)) && ((chr_range_high > VCF.pos(record) > chr_range_low))
                             push!(vcf_subarray,record)
                         end

                 else

                         if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos) && (VCF.hasfilter(record)) && (VCF.filter(record) == String["PASS"]) && ((VCF.chrom(record) == chr)) && ((chr_range_high > VCF.pos(record) > chr_range_low))
                             push!(vcf_subarray,record)

                         end
                 end
           end
    end

    return vcf_subarray
    end

"""
    pass_chrrange_filter(reader::GeneticVariation.VCF.Reader,sig_list,chr_range::AbstractString)
returns subarray of vcf records with io_pass_filter and io_chromosome_range_vcf_filter applied.
"""
    function pass_chrrange_filter(reader::GeneticVariation.VCF.Reader,chr_range::AbstractString)

        a=split(chr_range,":")
        chrwhole=a[1]
        chrnumber=split(chrwhole,"r")
        string_chr=chrnumber[2]
        chr=String(string_chr)
        range=a[2]
        splitrange=split(range, "-")
        lower_limit=splitrange[1]
        chr_range_low=parse(lower_limit)
        upper_limit=splitrange[2]
        chr_range_high=parse(upper_limit)

        vcf_subarray = Array{Any}(0)

               for record in reader

                      if typeof(VCF.chrom(record)) == String
                             chr = string(chr)

                             if (VCF.hasfilter(record)) && (VCF.filter(record) == String["PASS"]) && ((VCF.chrom(record) == chr)) && ((chr_range_high > VCF.pos(record) > chr_range_low))
                                 push!(vcf_subarray,record)
                             end

                     else

                             if (VCF.hasfilter(record)) && (VCF.filter(record) == String["PASS"]) && ((VCF.chrom(record) == chr)) && ((chr_range_high > VCF.pos(record) > chr_range_low))
                                 push!(vcf_subarray,record)

                             end
                     end
               end

        return vcf_subarray
    end

"""
        pass_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)
    returns subarray of vcf records with io_pass_filter, io_sig_list_vcf_filter, and io_chromosome_range_vcf_filter applied.
"""
    function pass_siglist_filter(vcf_filename,sig_list)

        vcf_subarray = Array{Any}(0)

        for row= 1:size(sig_list,1)
               dimension = size(sig_list,1)

               chr=(sig_list[row,1])
               pos=(sig_list[row,2])

               reader = VCF.Reader(open(vcf_filename, "r"))

               for record in reader

                      if typeof(VCF.chrom(record)) == String
                             chr = string(chr)

                             if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos) && (VCF.hasfilter(record)) && (VCF.filter(record) == String["PASS"])
                                 push!(vcf_subarray,record)
                             end

                     else

                             if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos) && (VCF.hasfilter(record)) && (VCF.filter(record) == String["PASS"])
                                 push!(vcf_subarray,record)

                             end
                     end
               end
        end

        return vcf_subarray
        end

"""
    chrrange_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)
returns subarray of vcf records with io_pass_filter, io_sig_list_vcf_filter, and io_chromosome_range_vcf_filter applied.
"""
        function chrrange_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)

            a=split(chr_range,":")
            chrwhole=a[1]
            chrnumber=split(chrwhole,"r")
            string_chr=chrnumber[2]
            chr=String(string_chr)
            range=a[2]
            splitrange=split(range, "-")
            lower_limit=splitrange[1]
            chr_range_low=parse(lower_limit)
            upper_limit=splitrange[2]
            chr_range_high=parse(upper_limit)

            vcf_subarray = Array{Any}(0)

            for row= 1:size(sig_list,1)
                   dimension = size(sig_list,1)

                   chr=(sig_list[row,1])
                   pos=(sig_list[row,2])

                   reader = VCF.Reader(open(vcf_filename, "r"))

                   for record in reader

                          if typeof(VCF.chrom(record)) == String
                                 chr = string(chr)

                                 if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos) && ((VCF.chrom(record) == chr)) && ((chr_range_high > VCF.pos(record) > chr_range_low))
                                     push!(vcf_subarray,record)
                                 end

                         else

                                 if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos) && ((VCF.chrom(record) == chr)) && ((chr_range_high > VCF.pos(record) > chr_range_low))
                                     push!(vcf_subarray,record)

                                 end
                         end
                   end
            end

            return vcf_subarray
            end

#functions for converting vcf record array to numerical array

"""
    create_chr_dict()
creates dict for use in combined_all_genotype_array_functions() for removing 'chr' from chromosome labels to allow sorting variant records by chromosome position.
"""
function create_chr_dict()

    chr_dict = Dict()

    labels_with_chr = ["chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM"]

    for item in labels_with_chr
           chr_dict[item] = 800
    end

    return chr_dict
end

    #=
    replace_chr_dict = Dict(r"(chr)" => "")
    join([replace_chr_dict[c] for c in genotype_array[:,1]])
    =#

"""
    combined_all_genotype_array_functions(sub)
convert sub from variant filters to gt_num_array and gt_chromosome_labels for plot functions.
"""
function combined_all_genotype_array_functions(sub)
    genotype_array = generate_genotype_array(sub,"GT")

    map!(s->replace(s, "chr", ""), genotype_array, genotype_array)

    #insert remove_chr_col1_of_genotype_array - no need to clean siglist - only support human v0.1
    #insert function to replace chr here

#=

chr_dict = create_chr_dict()

function translate(c)
       chr_dict[c]
end

no_chr_genotype_array = map(translate, genotype_array)

    if ismatch(r"^chr",genotype_array[1])
        println("replacing chr in chr_labels")

        for i in genotype_array[:,1]
            println(i)
            n=replace(i, "chr" => "")
            println("$n")
        end
    end
=#

    clean_column1!(genotype_array)
    genotype_array=ViVa.sort_genotype_array(genotype_array)
    geno_dict = define_geno_dict()
    gt_num_array,gt_chromosome_labels = translate_genotype_to_num_array(genotype_array, geno_dict)

    return gt_num_array,gt_chromosome_labels
end

"""
    combined_all_read_depth_array_functions(sub)
convert sub from variant filters to dp_num_array and dp_chromosome_labels for plot functions.
"""
function combined_all_read_depth_array_functions(sub)

    read_depth_array = ViVa.generate_genotype_array(sub,"DP")
    map!(s->replace(s, "chr", ""), read_depth_array, read_depth_array)
    clean_column1!(read_depth_array)
    read_depth_array=ViVa.sort_genotype_array(read_depth_array)
    dp_num_array,dp_chromosome_labels = translate_readdepth_strings_to_num_array(read_depth_array)

    return dp_num_array,dp_chromosome_labels
end

"""
    generate_genotype_array(record_sub::Array{Any,1},genotype_field::String)
Returns numerical array of genotype values (either genotype or read_depth values) which are translated by another function into num_array
Where genotype_field is either GT or DP to visualize genotype or read_depth
"""
function generate_genotype_array(record_sub::Array{Any,1},y)

       num_samples = length(VCF.genotype(record_sub[1]))

       df = DataFrame(String, 0, num_samples+2)

       for record in record_sub

              genotype_data_per_variant = VCF.genotype(record, 1:num_samples, y)
              vcf_chrom = VCF.chrom(record)
              vcf_pos = string(VCF.pos(record))
              genotype_data_per_variant = vcat(vcf_pos,genotype_data_per_variant)
              genotype_data_per_variant = vcat(vcf_chrom,genotype_data_per_variant)
              push!(df, genotype_data_per_variant)
       end

       pre_num_array = Matrix(df)

       return pre_num_array
end

"""
    define_geno_dict()
returns dictionary of values for use in replace_genotype_with_vals()
"""
function define_geno_dict()
    geno_dict = Dict()

    homo_variant = ["1/1" "1/2" "2/2" "1/3" "2/3" "3/3" "1/4" "2/4" "3/4" "4/4" "1/5" "2/5" "3/5" "4/5" "5/5" "1/6" "2/6" "3/6" "4/6" "5/6" "6/6" "1|1" "1|2" "2|2" "1|3" "2|3" "3|3" "1|4" "2|4" "3|4" "4|4" "1|5" "2|5" "3|5" "4|5" "5|5" "1|6" "2|6" "3|6" "4|6" "5|6" "6|6" "2/1" "3/2" "4/2" "5/2" "6/2" "4/3" "5/3" "6/3" "5/4" "6/4" "6/5" "2|1" "3|2" "4|2" "5|2" "6|2" "4|3" "5|3" "6|3" "5|4" "6|4" "6|5"]

    hetero_variant = ["0/1" "0/2" "0/3" "0/4" "0/5" "0/6" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "0|1" "0|2" "0|3" "0|4" "0|5" "0|6" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0"]

    no_data = ["./." ".|."]

    ref = ["0/0" "0|0"]

    for item in homo_variant
           geno_dict[item] = 800
    end
    for item in hetero_variant
           geno_dict[item] = 600
    end
    for item in no_data
           geno_dict[item] = 0
    end
    for item in ref
           geno_dict[item] = 400
    end

    return geno_dict
end

"""
    translate_genotype_to_num_array(genotype_array,geno_dict)
returns a tuple of num_array for plotting, and chromosome labels for plotting as label bar.
Translates array of genotype values to numerical array of categorical values.
Genotype values are converted to categorical values. No_call=0, 0/0=400, heterozygous_variant=600, homozygous_variant=800
"""
function translate_genotype_to_num_array(genotype_array,geno_dict)

    array_for_plotly = genotype_array[:,3:size(genotype_array,2)]
    chromosome_labels = genotype_array[:,1:2]

        function translate(c)
               geno_dict[c]
        end

    num_array = map(translate, array_for_plotly)

    return num_array,chromosome_labels
end

"""
    translate_readdepth_strings_to_num_array(read_depth_array::Array{Any,2})
Returns array of read_depth as int for plotting and average calculation.
By default, read depth values over 100 are replaced with 100 to improve heatmap visualization (see read_depth_threshhold() ).
Where read_depth_array is output of generate_genotype_array() for DP option
returns a tuple of num_array type Int for average calculation and plotting, and chromosome labels for plotting as label bar
"""
function translate_readdepth_strings_to_num_array(read_depth_array::Array{Any,2})

       chromosome_labels = read_depth_array[:,1:2]

       dp_array_for_plotly = read_depth_array[:,3:size(read_depth_array,2)]

       map!(s->replace(s, ".", "0"), dp_array_for_plotly, dp_array_for_plotly)

       dp_array_for_plotly = [parse(Int, i) for i in dp_array_for_plotly]

       for i in dp_array_for_plotly
              if i == "."
                     #println(".")
              end
       end

       return dp_array_for_plotly, chromosome_labels
end


#functions for sample filters

"""
    get_sample_names(reader)
returns sample ids of vcf file as a vector of symbols for naming columns of num_array dataframe object for column filter functions
"""
function get_sample_names(reader)

    hdr=VCF.header(reader)
    sample_names=hdr.sampleID
    sample_names=reshape(sample_names,1,length(sample_names))
    sample_names=convert(Array{Symbol}, sample_names)

    return sample_names
end

"""
    find_group_label_indices(pheno)
find indices and determines names for group 1 and group 2 labels on plots. finds index of center of each sample group to place tick mark and label.
"""
function find_group_label_indices(pheno,trait_to_group_by,row_to_sort_by)

    group1_label=split(trait_to_group_by,",")[1]
    group2_label=split(trait_to_group_by,",")[2]
    group1_index = ((findlast(pheno[row_to_sort_by,:],1)) - (findfirst(pheno[row_to_sort_by,:],1)))/2
    group2_index = ((findlast(pheno[row_to_sort_by,:],2) - (findfirst(pheno[row_to_sort_by,:],2)))/2) + (findlast(pheno[row_to_sort_by,:],1)) + 0.5
    group_dividing_line = findfirst(pheno[row_to_sort_by,:],2) - (0.5)

    group_label_pack = [group1_index, group2_index, group_dividing_line, group1_label, group2_label]

    return group_label_pack
end

"""
    sortcols_by_phenotype_matrix(pheno_matrix_filename::String,trait_to_group_by::String,num_array::Array{Int64,2}, sample_names::Array{Symbol,2})
group samples by a common trait using a user generated key matrix ("phenotype matrix")
"""
function sortcols_by_phenotype_matrix(pheno_matrix_filename::String,trait_to_group_by::String, num_array::Array{Int64,2}, sample_names::Array{Symbol,2})

    pheno = readdlm(pheno_matrix_filename, ',')

    #get row numer of user-chosen phenotype characteristic to sort columns by
    row_to_sort_by = find(x -> x == trait_to_group_by, pheno)
    row_to_sort_by = row_to_sort_by[1]

    #remove phenotype_row_labels used to identify row to sort by, so row can be sorted without strings causing issues
    pheno = pheno[:,2:size(pheno,2)]

    pheno = sortcols(pheno, by = x -> x[row_to_sort_by], rev = false)

    id_list = pheno[1,:]

    #vcf_header = names(df_vcf)
    #vcf_info_columns = vcf_header[1:9]

    sample_ids=sample_names

    col_new_order=vec(id_list)

    col_new_order = [Symbol(i) for i in col_new_order]

    df1_vcf = DataFrame(num_array)

    rename!(df1_vcf, f => t for (f, t) = zip(names(df1_vcf), sample_ids))

    #println(typeof(names(df1_vcf))) #Array{Symbol,1} | same in test
    #println(typeof(col_new_order)) #Array{Any,1} | same in test

    vcf = df1_vcf[:, col_new_order]

    vcf = Matrix(vcf)

    group_label_pack = find_group_label_indices(pheno,trait_to_group_by,row_to_sort_by)

    return vcf,group_label_pack
end

"""
    select_columns(filename_sample_list::AbstractString, num_array::Array{Int64,2}, sample_names::Array{Symbol,2})
returns num_array with columns matching user generated list of sample ids to select for analysis. num_array now has sample ids in first row.
"""
function select_columns(filename_sample_list::AbstractString, num_array::Array{Int64,2}, sample_names::Array{Symbol,2})

    selectedcolumns=readdlm(filename_sample_list)

    header_as_strings = selectedcolumns

    for item = 1:size(selectedcolumns,2)
        #selectedcolumns[item] = replace(selectedcolumns[item],".","_")  #names in sample list dont match vcf so have to clean
        selectedcolumns[item] = Symbol(selectedcolumns[item])
    end

    col_selectedcolumns=vec(selectedcolumns)

    df_num_array = DataFrame(num_array)

    rename!(df_num_array, f => t for (f, t) = zip(names(df_num_array), sample_names))

    df_selected_samples_num_array = df_num_array[:, col_selectedcolumns]

    selected_samples_num_array = Matrix(df_selected_samples_num_array)
    println(filename_sample_list)
    return selected_samples_num_array
end


#functions for mathematic analysis
"""
    avg_dp_samples(dp_num_array::Array{Int64,2})
create sample_avg_list vector that lists averages of read depth for each sample for input into avg_sample_dp_line_chart(sample_avg_list)
dp_num_array must contain dp values as Int64 and be without chromosome position columns
"""
function avg_dp_samples(dp_num_array::Array{Int64,2})

    avg_dps_all = Array{Float64}(0)

    for column = 1:size(dp_num_array,2)

        all_dps_patient = dp_num_array[1:size(dp_num_array,1),column]
        sum_dp = sum(all_dps_patient)
        avg_dp = sum_dp/(size(dp_num_array,1))
        push!(avg_dps_all,avg_dp)
    end

    return avg_dps_all
end


"""
    avg_dp_variant(dp_num_array::Array{Int64,2})
create variant_avg_list vector that lists averages of read depth for each variant for input into avg_variant_dp_line_chart(variant_avg_list)
"""
function avg_dp_variant(dp_num_array::Array{Int64,2})

    avg_dps_all = Array{Float64}(0)

    for row = 1:size(dp_num_array,1)

        all_dps_variant = dp_num_array[row,1:size(dp_num_array,2)]
        sum_dp = sum(all_dps_variant)
        avg_dp = sum_dp/(size(dp_num_array,2))
        push!(avg_dps_all,avg_dp)
    end

    return avg_dps_all
end

"""
    list_sample_names_low_dp(sample_avg_list::Array{Float64,2},sample_names)
returns list of sample ids that have an average read depth of under 15 across all variant positions
"""
function list_sample_names_low_dp(sample_avg_list::Array{Float64,1},sample_names)

    low_dp_index_list = Array{Int64}(0)

        for item = 1:length(sample_avg_list)
            if sample_avg_list[item] < 15
                push!(low_dp_index_list,item)
            end
        end

    low_dp_sample_names = Array{String}(0)

        for i in low_dp_index_list
            sample_name = sample_names[i]
            string_sample_name=string(sample_name)
            push!(low_dp_sample_names,string_sample_name)
        end

        return low_dp_sample_names
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

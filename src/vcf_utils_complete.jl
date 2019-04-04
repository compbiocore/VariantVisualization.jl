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
        matrix_with_chr_column[i, 1] = matrix_with_chr_column[i, 1] == "MT" ? "25" : matrix_with_chr_column[i, 1]
        matrix_with_chr_column[i, 1] = matrix_with_chr_column[i, 1] == "Mt" ? "25" : matrix_with_chr_column[i, 1]
    end
end

"""
    clean_column1_siglist!(siglist)
Replaces strings "X","Y","M" with 23,24,25 {Int} in array generated in load_siglist()
use in load_siglist() because X and Y need to be replaced with Int
"""
function clean_column1_siglist!(siglist)

    for i = 1:(size(siglist, 1))
        siglist[i, 1] = siglist[i, 1] == "X" ? 23 : siglist[i, 1]
        siglist[i, 1] = siglist[i, 1] == "Y" ? 24 : siglist[i, 1]
        siglist[i, 1] = siglist[i, 1] == "M" ? 25 : siglist[i, 1]
        siglist[i, 1] = siglist[i, 1] == "MT" ? 25 : siglist[i, 1]
        siglist[i, 1] = siglist[i, 1] == "Mt" ? 25 : siglist[i, 1]
    end
end

"""
    returnXY_column1!(chr_label_vector)
Replace String "23","24","25" with "X","Y","M" in chromosome label vector used for plot labels
"""
function returnXY_column1!(chr_label_vector)

    for i = 1:(size(chr_label_vector, 1))
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == 23 ? "X" : chr_label_vector[i, 1]
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == 24 ? "Y" : chr_label_vector[i, 1]
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == 25 ? "M" : chr_label_vector[i, 1]
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == 25 ? "MT" : chr_label_vector[i, 1]
        chr_label_vector[i, 1] = chr_label_vector[i, 1] == 25 ? "Mt" : chr_label_vector[i, 1]
    end
end

"""
    returnXY_column1_siglist!(siglist_sorted)
Replace String "23","24","25" with "X","Y","M" in siglist for filtering
"""
function returnXY_column1_siglist!(siglist_sorted)

    for i = 1:(size(siglist_sorted, 1))
        siglist_sorted[i, 1] = siglist_sorted[i, 1] == 23 ? "X" : siglist_sorted[i, 1]
        siglist_sorted[i, 1] = siglist_sorted[i, 1] == 24 ? "Y" : siglist_sorted[i, 1]
        siglist_sorted[i, 1] = siglist_sorted[i, 1] == 25 ? "M" : siglist_sorted[i, 1]
        siglist_sorted[i, 1] = siglist_sorted[i, 1] == 25 ? "MT" : siglist_sorted[i, 1]
        siglist_sorted[i, 1] = siglist_sorted[i, 1] == 25 ? "Mt" : siglist_sorted[i, 1]
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
    genotype_array = sortslices(genotype_array, dims=1, by=x->(x[1],x[2]))

return genotype_array
end

"""
    load_siglist(filename::AbstractString)
where x = filename of significant SNP variant location list in comma delimited format (saved as .csv)
"""
function load_siglist(filename::AbstractString)

siglist_unsorted = readdlm(filename, ',', skipstart=1)

VariantVisualization.clean_column1_siglist!(siglist_unsorted)
#siglist = sortrows(siglist_unsorted, by=x->(x[1],x[2]))#version
siglist = sortslices(siglist_unsorted, dims=1, by=x->(x[1],x[2]))

returnXY_column1_siglist!(siglist)

if typeof(siglist) == Array{Float64,2}
    siglist = trunc.(Int,siglist)
    return siglist
end

return siglist

end

"""
    build_set_from_list(sig_list::Array{Any,2})
build set of tuples of chrom and pos of each record in vcf for use in sig_list_filters.
"""
function build_set_from_list(sig_list::Array{Any,2})
    position_set = Set{Tuple{Any,Int}}()  # For now, assumes all chr names are Int...

    for row in 1:size(sig_list,1)
        chr = sig_list[row,1]
        pos = sig_list[row,2]
        push!(position_set, (string(chr), Int(pos)))
    end

    return position_set
end

#functions for variant filters
"""
io_genomic_range_vcf_filter(chr_range::String, vcf_filename::AbstractString)
create subarray of vcf variant records matching user specified chromosome range in format: (e.g. chr1:0-30000000)
"""
function io_genomic_range_vcf_filter(chr_range::String,vcf_filename::AbstractString)
    a=split(chr_range,":")
    chrwhole=a[1]
    chrnumber=split(chrwhole,"r")
    string_chr=chrnumber[2]
    chr=String(string_chr)
    range=a[2]
    splitrange=split(range, "-")
    lower_limit=splitrange[1]
    chr_range_low=Meta.parse(lower_limit)
    upper_limit=splitrange[2]
    chr_range_high=Meta.parse(upper_limit)

    reader = VCF.Reader(open(vcf_filename, "r"))
    record1=VCF.read(reader)

    open(VCF.Reader, vcf_filename) do reader

    vcf_record = VCF.Record()
    vcf_subarray = Array{Any}(undef,0)

    if occursin("chr",VCF.chrom(record1))

        while !eof(reader)
            read!(reader, vcf_record)

               if (VCF.chrom(vcf_record) == chrwhole) && (chr_range_high >= VCF.pos(vcf_record) >= chr_range_low)
                      push!(vcf_subarray,copy(vcf_record))
               end
        end

    else
        while !eof(reader)
            read!(reader, vcf_record)

               if (VCF.chrom(vcf_record) == chr) && (chr_range_high >= VCF.pos(vcf_record) >= chr_range_low)
                      push!(vcf_subarray,copy(vcf_record))
               end
        end
    end
    return vcf_subarray
end
end

#=
function io_genomic_range_vcf_filter(chr_range::String,vcf_filename::AbstractString)

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

       vcf_subarray = Array{Any}(undef,0)

       reader = VCF.Reader(open(vcf_filename, "r"))
       record1=VCF.read(reader)

       #if occursin("chr",VCF.chrom(record1))#version
       if occursin("chr",VCF.chrom(record1))
           for record in reader

                  if (VCF.chrom(record) == chrwhole) && (chr_range_high > VCF.pos(record) > chr_range_low)
                         push!(vcf_subarray,record)
                  end
           end

       else
           for record in reader

                  if (VCF.chrom(record) == chr) && (chr_range_high > VCF.pos(record) > chr_range_low)
                         push!(vcf_subarray,record)
                  end
           end
       end

       return vcf_subarray
end
=#

"""
    io_sig_list_vcf_filter(sig_list,vcf_filename)
returns subarray of variant records matching a list of variant positions returned from load_siglist()
"""
function io_sig_list_vcf_filter(sig_list, vcf_filename)
        open(VCF.Reader, vcf_filename) do reader

        vcf_subarray = Array{Any}(undef,0)

        posset = build_set_from_list(sig_list)
        vcf_record = VCF.Record()

        while !eof(reader)
            read!(reader, vcf_record)
            chr = VCF.chrom(vcf_record)
            pos = VCF.pos(vcf_record)

            if (chr, pos) in posset
                push!(vcf_subarray, copy(vcf_record))
            end
        end

        return vcf_subarray
    end
end

"""
    io_pass_filter(vcf_filename)
returns subarray of vcf records including only records with FILTER status = PASS
"""
function io_pass_filter(vcf_filename)

    open(VCF.Reader, vcf_filename) do reader

    vcf_subarray = Array{Any}(undef,0)

    vcf_record = VCF.Record()

    while !eof(reader)
        read!(reader, vcf_record)

       if VCF.hasfilter(vcf_record) && VCF.filter(vcf_record) == String["PASS"]
              push!(vcf_subarray,copy(vcf_record))
       end
end

return vcf_subarray
end
end

"""
    pass_genomic_range_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)
returns subarray of vcf records with io_pass_filter, io_sig_list_vcf_filter, and io_genomic_range_vcf_filter applied.
"""
function pass_genomic_range_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)

    a=split(chr_range,":")
    chrwhole=a[1]
    chrnumber=split(chrwhole,"r")
    string_chr=chrnumber[2]
    chr=String(string_chr)
    range=a[2]
    splitrange=split(range, "-")
    lower_limit=splitrange[1]
    chr_range_low=Meta.parse(lower_limit)
    upper_limit=splitrange[2]
    chr_range_high=Meta.parse(upper_limit)

    posset = build_set_from_list(sig_list)

    vcf_subarray = Array{Any}(undef,0)
    vcf_record = VCF.Record()

    reader = VCF.Reader(open(vcf_filename, "r"))
    record1=VCF.read(reader)

    open(VCF.Reader, vcf_filename) do reader

    while !eof(reader)
        read!(reader, vcf_record)

        if size(vcf_subarray,1) == size(sig_list,1)
                       break
        end

    for row = 1:size(sig_list,1)
           dimension = size(sig_list,1)

           chr_sig=(sig_list[row,1])
           pos_sig=(sig_list[row,2])

                  if typeof(VCF.chrom(vcf_record)) == String

                      chr_sig = string(chr_sig)

                         #if occursin("chr",VCF.chrom(record1))#version
                         if occursin("chr",VCF.chrom(record1))

                             if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chrwhole)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                 push!(vcf_subarray,copy(vcf_record))
                                 break
                             end

                        else

                            if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chr)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                push!(vcf_subarray,copy(vcf_record))
                                break
                            end
                        end

                 else

                     if occursin("chr",VCF.chrom(record1))

                         if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chrwhole)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                             push!(vcf_subarray,copy(vcf_record))
                             break
                         end

                     else
                        if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chr)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                            push!(vcf_subarray,copy(vcf_record))
                            break
                        end
                     end

                 end
             end
           end
    end

    return vcf_subarray
    end

"""
    pass_genomic_range_filter(reader::GeneticVariation.VCF.Reader,chr_range::AbstractString,vcf_filename)
returns subarray of vcf records with io_pass_filter and io_genomic_range_vcf_filter applied.
"""
    function pass_genomic_range_filter(reader::GeneticVariation.VCF.Reader,chr_range::AbstractString,vcf_filename)

        a=split(chr_range,":")
        chrwhole=a[1]
        chrnumber=split(chrwhole,"r")
        string_chr=chrnumber[2]
        chr=String(string_chr)
        range=a[2]
        splitrange=split(range, "-")
        lower_limit=splitrange[1]
        chr_range_low=Meta.parse(lower_limit)
        upper_limit=splitrange[2]
        chr_range_high=Meta.parse(upper_limit)

        vcf_subarray = Array{Any}(undef,0)

        reader = VCF.Reader(open(vcf_filename, "r"))
        record1=VCF.read(reader)

        vcf_record = VCF.Record()

        open(VCF.Reader, vcf_filename) do reader

            while !eof(reader)
                read!(reader, vcf_record)

                 if occursin("chr",VCF.chrom(record1))

                      if typeof(VCF.chrom(vcf_record)) == String

                             chr = string(chr)

                             if (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chrwhole)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                 push!(vcf_subarray,copy(vcf_record))

                             end

                     else

                             if (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chrwhole)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                 push!(vcf_subarray,copy(vcf_record))
                             end
                     end


                 else

                     if typeof(VCF.chrom(vcf_record)) == String
                            chr = string(chr)

                            if (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chr)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                push!(vcf_subarray,copy(vcf_record))
                            end

                    else

                            if (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"]) && ((VCF.chrom(vcf_record) == chr)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                push!(vcf_subarray,copy(vcf_record))

                            end
                    end

               end
           end
       end

        return vcf_subarray
    end

"""
    pass_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)
returns subarray of vcf records with io_pass_filter, io_sig_list_vcf_filter, and io_genomic_range_vcf_filter applied.
"""
function pass_siglist_filter(vcf_filename,sig_list)

    vcf_subarray = Array{Any}(undef,0)
    open(VCF.Reader, vcf_filename) do reader

    posset = build_set_from_list(sig_list)
    vcf_record = VCF.Record()

    while !eof(reader)
        read!(reader, vcf_record)

        if size(vcf_subarray,1) == size(sig_list,1)
                       break
        end

        for row= 1:size(sig_list,1)
               dimension = size(sig_list,1)

               chr=(sig_list[row,1])
               pos=(sig_list[row,2])

                      if typeof(VCF.chrom(vcf_record)) == String
                             chr = string(chr)

                             if (VCF.chrom(vcf_record) == chr) && (VCF.pos(vcf_record) == pos) && (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"])
                                push!(vcf_subarray, copy(vcf_record))
                                 break
                             end

                     else

                             if (VCF.chrom(vcf_record) == chr) && (VCF.pos(vcf_record) == pos) && (VCF.hasfilter(vcf_record)) && (VCF.filter(vcf_record) == String["PASS"])
                                push!(vcf_subarray, copy(vcf_record))
                                 break
                             end
                     end
        end
    end
end

    return vcf_subarray
end

"""
    genomic_range_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)
returns subarray of vcf records with io_pass_filter, io_sig_list_vcf_filter, and io_genomic_range_vcf_filter applied.
"""
function genomic_range_siglist_filter(vcf_filename,sig_list,chr_range::AbstractString)

    a=split(chr_range,":")
    chrwhole=a[1]
    chrnumber=split(chrwhole,"r")
    string_chr=chrnumber[2]
    chr=String(string_chr)
    range=a[2]
    splitrange=split(range, "-")
    lower_limit=splitrange[1]
    chr_range_low=Meta.parse(lower_limit)
    upper_limit=splitrange[2]
    chr_range_high=Meta.parse(upper_limit)

    vcf_subarray = Array{Any}(undef,0)

    reader = VCF.Reader(open(vcf_filename, "r"))
    record1=VCF.read(reader)

    open(VCF.Reader, vcf_filename) do reader

    vcf_record = VCF.Record()

    while !eof(reader)
        read!(reader, vcf_record)

        if size(vcf_subarray,1) == size(sig_list,1)
                       break
        end

        for row= 1:size(sig_list,1)
            dimension = size(sig_list,1)

            chr_sig=(sig_list[row,1])
            pos_sig=(sig_list[row,2])

            if occursin("chr",VCF.chrom(record1))

                    if typeof(VCF.chrom(vcf_record)) == String
                             chr_sig = string(chr_sig)

                             if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && ((VCF.chrom(vcf_record) == chrwhole)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                 push!(vcf_subarray,copy(vcf_record))
                                 break
                             end

                     else

                             if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && ((VCF.chrom(vcf_record) == chrwhole)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                 push!(vcf_subarray,copy(vcf_record))
                                 break
                             end
                     end

            else

                     if typeof(VCF.chrom(vcf_record)) == String
                              chr_sig = string(chr_sig)

                              if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && ((VCF.chrom(vcf_record) == chr)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                  push!(vcf_subarray,copy(vcf_record))
                                  break
                              end

                      else

                              if (VCF.chrom(vcf_record) == chr_sig) && (VCF.pos(vcf_record) == pos_sig) && ((VCF.chrom(vcf_record) == chr)) && ((chr_range_high >= VCF.pos(vcf_record) >= chr_range_low))
                                  push!(vcf_subarray,copy(vcf_record))
                                  break
                              end
                      end
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

    #labels_with_chr = []

    chr_dict = Dict("chr1"=> 1,"chr2"=> 2,"chr3"=> 3,"chr4"=> 4,"chr5"=> 5,"chr6"=> 6,"chr7"=> 7,"chr8"=> 8,"chr9"=> 9,"chr10"=> 10,"chr11"=> 11,"chr12"=> 12,"chr13"=> 13,"chr14"=> 14,"chr15"=> 15,"chr16"=> 16,"chr17"=> 17,"chr18"=> 18,"chr19"=> 19,"chr20"=> 20,"chr21"=> 21,"chr22"=> 22,"chrX"=> X,"chrY"=> Y,"chrM"=> M)

    return chr_dict
end

"""
    combined_all_genotype_array_functions(sub)
convert sub from variant filters to gt_num_array and gt_chromosome_labels for plot functions.
"""
function combined_all_genotype_array_functions(sub)
    genotype_array = generate_genotype_array(sub,"GT")
    map!(s->replace(s, "chr" => ""), genotype_array, genotype_array)
    clean_column1!(genotype_array)
    genotype_array=VariantVisualization.sort_genotype_array(genotype_array)
    geno_dict = define_geno_dict()
    gt_num_array,gt_chromosome_labels = translate_genotype_to_num_array(genotype_array, geno_dict)

    return gt_num_array,gt_chromosome_labels
end

"""
    combined_all_read_depth_array_functions(sub)
convert sub from variant filters to dp_num_array and dp_chromosome_labels for plot functions.
"""
function combined_all_read_depth_array_functions(sub)

    read_depth_array = VariantVisualization.generate_genotype_array(sub,"DP")
    map!(s->replace(s, "chr"=>""), read_depth_array, read_depth_array)
    clean_column1!(read_depth_array)
    read_depth_array=VariantVisualization.sort_genotype_array(read_depth_array)
    dp_num_array,dp_chromosome_labels = translate_readdepth_strings_to_num_array(read_depth_array)
    return dp_num_array,dp_chromosome_labels
end

"""
    combined_all_read_depth_array_functions_for_avg_dp(sub)
convert sub from variant filters to dp_num_array and dp_chromosome_labels for plot functions.
"""
function combined_all_read_depth_array_functions_for_avg_dp(sub)

    read_depth_array = VariantVisualization.generate_genotype_array(sub,"DP")
    map!(s->replace(s, "chr"=>""), read_depth_array, read_depth_array)
    clean_column1!(read_depth_array)
    read_depth_array=VariantVisualization.sort_genotype_array(read_depth_array)
    dp_num_array,dp_chromosome_labels = translate_readdepth_strings_to_num_array_for_avg_dp(read_depth_array)
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

    homo_variant = ["1/1" "1/2" "2/2" "1/3" "2/3" "3/3" "1/4" "2/4" "3/4" "4/4" "1/5" "2/5" "3/5" "4/5" "5/5" "1/6" "2/6" "3/6" "4/6" "5/6" "6/6" "1|1" "1|2" "2|2" "1|3" "2|3" "3|3" "1|4" "2|4" "3|4" "4|4" "1|5" "2|5" "3|5" "4|5" "5|5" "1|6" "2|6" "3|6" "4|6" "5|6" "6|6" "2/1" "3/1" "4/1" "5/1" "6/1" "3/2" "4/2" "5/2" "6/2" "4/3" "5/3" "6/3" "5/4" "6/4" "6/5" "2|1" "3|1" "4|1" "5|1" "6|1" "3|2" "4|2" "5|2" "6|2" "4|3" "5|3" "6|3" "5|4" "6|4" "6|5"]

    hetero_variant = ["0/1" "0/2" "0/3" "0/4" "0/5" "0/6" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "0|1" "0|2" "0|3" "0|4" "0|5" "0|6" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0"]

    no_data = ["./." ".|."]

    ref = ["0/0" "0|0"]

    for item in homo_variant
           geno_dict[item] = 3
    end
    for item in hetero_variant
           geno_dict[item] = 2
    end
    for item in no_data
           geno_dict[item] = 0
    end
    for item in ref
           geno_dict[item] = 1
    end

    return geno_dict
end

"""
    translate_genotype_to_num_array(genotype_array,geno_dict)
returns a tuple of num_array for plotting, and chromosome labels for plotting as label bar.
Translates array of genotype values to numerical array of categorical values.
Genotype values are converted to categorical values. No_call=0, 0/0=1, heterozygous_variant=2, homozygous_variant=3
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
Ceiling of dp=100 is set to prevent high dp value from hiding (or "blowing out") low dp values. (see read_depth_threshhold() ).
Where read_depth_array is output of generate_genotype_array() for DP option
returns a tuple of num_array type Int for average calculation and plotting, and chromosome labels for plotting as label bar
"""
function translate_readdepth_strings_to_num_array(read_depth_array::Array{Any,2})

       chromosome_labels = read_depth_array[:,1:2]

       dp_array_for_plotly = read_depth_array[:,3:size(read_depth_array,2)]

       map!(s->replace(s, "."=>"-20"), dp_array_for_plotly, dp_array_for_plotly)

       dp_array_for_plotly = [parse(Int, i) for i in dp_array_for_plotly]

       return dp_array_for_plotly, chromosome_labels
end

"""
    translate_readdepth_strings_to_num_array_for_avg_dp(read_depth_array::Array{Any,2})
Returns array of read_depth as int for plotting and average calculation.
'read_depth_array' is output of generate_genotype_array() for DP option
returns a tuple of num_array type Int for average calculation and plotting, and chromosome labels for plotting as label bar
No call is replaced with 0 for avg_calculation. Ceiling of dp=100 is set to prevent high dp value from hiding (or "blowing out") low dp values.
"""
function translate_readdepth_strings_to_num_array_for_avg_dp(read_depth_array::Array{Any,2})

       chromosome_labels = read_depth_array[:,1:2]

       dp_array_for_plotly = read_depth_array[:,3:size(read_depth_array,2)]

       map!(s->replace(s, "."=>"0"), dp_array_for_plotly, dp_array_for_plotly)

       dp_array_for_plotly = [parse(Int, i) for i in dp_array_for_plotly]

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

    #sample_names=convert(Array{Symbol}, sample_names)#version
    sample_names = [Symbol(i) for i in sample_names]

    return sample_names
end

"""
    find_group_label_indices(pheno,trait_to_group_by,row_to_sort_by)
find indices and determines names for group 1 and group 2 labels on plots. finds index of center of each sample group to place tick mark and label.
"""
function find_group_label_indices(pheno,trait_to_group_by,row_to_sort_by)

    group1_label=split(trait_to_group_by,",")[1]
    group2_label=split(trait_to_group_by,",")[2]

    group1_index = (something(findlast(isequal(1),pheno[row_to_sort_by,:]),0) -
    something(findfirst(isequal(1),pheno[row_to_sort_by,:]),0))/2

    group2_index = (something(findlast(isequal(2),pheno[row_to_sort_by,:]),0) -
    something(findfirst(isequal(2),pheno[row_to_sort_by,:]),0))/2 +
    something(findlast(isequal(1),pheno[row_to_sort_by,:]),0) + 0.5

    group_dividing_line = something(findfirst(isequal(2),pheno[row_to_sort_by,:]),0) - 0.5

    #group1_index = ((findlast(pheno[row_to_sort_by,:],1)) - (findfirst(pheno[row_to_sort_by,:],1)))/2 #version
    #group2_index = ((findlast(pheno[row_to_sort_by,:],2) - (findfirst(pheno[row_to_sort_by,:],2)))/2) + (findlast(pheno[row_to_sort_by,:],1)) + 0.5#version
    #group_dividing_line = findfirst(pheno[row_to_sort_by,:],2) - (0.5)#version

    group_label_pack = [group1_index, group2_index, group_dividing_line, group1_label, group2_label]

    return group_label_pack
end

"""
    sortcols_by_phenotype_matrix(pheno_matrix_filename::String,trait_to_group_by::String,num_array::Array{Int64,2}, sample_names::Array{Symbol,2})
group samples by a common trait using a user generated key matrix ("phenotype matrix") returns num_array,group_label_pack,
"""
function sortcols_by_phenotype_matrix(pheno_matrix_filename::String,trait_to_group_by::String, num_array::Array{Int64,2}, sample_names)

    pheno = readdlm(pheno_matrix_filename, ',')

    #get row numer of user-chosen phenotype characteristic to sort columns by
    row_to_sort_by = findall(x -> x == trait_to_group_by, pheno)
    row_to_sort_by = row_to_sort_by[1]

    #remove phenotype_row_labels used to identify row to sort by, so row can be sorted without strings causing issues
    pheno_no_trait_labels = pheno[:,2:size(pheno,2)]
    trait_labels=pheno[2:size(pheno,1),1]

    #pheno_no_trait_labels = sortcols(pheno_no_trait_labels, by = x -> x[row_to_sort_by], rev = false)#version
    pheno_no_trait_labels = sortslices(pheno_no_trait_labels, dims=2, by = x -> x[row_to_sort_by], rev = false)

    id_list = pheno_no_trait_labels[1,:]
    #pheno_data=pheno_no_trait_labels[]#reshape

    sample_ids=vec(sample_names)

    col_new_order=vec(id_list)

    col_new_order = [Symbol(i) for i in col_new_order]

    df1_vcf = DataFrame(num_array)

    #rename!(df1_vcf, f => t for (f, t) = zip(names(df1_vcf), sample_ids))#version
    names!(df1_vcf, sample_ids)

    ordered_num_array = df1_vcf[:, col_new_order]

    ordered_num_array = Matrix(ordered_num_array)

    group_label_pack = find_group_label_indices(pheno_no_trait_labels,trait_to_group_by,row_to_sort_by)

    id_list=reshape(id_list,1,length(id_list))

    return ordered_num_array,group_label_pack,pheno_no_trait_labels,id_list,trait_labels

end

"""
    select_columns(filename_sample_list::AbstractString, num_array::Array{Int64,2}, sample_names)
returns num_array with columns matching user generated list of sample ids to select for analysis. num_array now has sample ids in first row.
"""
function select_columns(filename_sample_list::AbstractString, num_array::Array{Int64,2}, sample_names)

    selectedcolumns=readdlm(filename_sample_list)

    header_as_strings = selectedcolumns

    for item = 1:size(selectedcolumns,2)
        selectedcolumns[item] = Symbol(selectedcolumns[item])
    end

    col_selectedcolumns=vec(selectedcolumns)

    #selectedcolumns=convert(Array{Symbol}, col_selectedcolumns) #version
    selectedcolumns=Symbol(col_selectedcolumns)

    df_num_array = DataFrame(num_array)

    rename!(df_num_array, f => t for (f, t) = zip(names(df_num_array), sample_names))

    df_selected_samples_num_array = df_num_array[:, col_selectedcolumns]

    selected_samples_num_array = Matrix(df_selected_samples_num_array)

    return selected_samples_num_array,col_selectedcolumns
end


#functions for mathematic analysis
"""
    avg_dp_samples(dp_num_array::Array{Int64,2})
create sample_avg_list vector that lists averages of read depth for each sample for input into avg_sample_dp_line_chart(sample_avg_list)
dp_num_array must contain dp values as Int64 and be without chromosome position columns
"""
function avg_dp_samples(dp_num_array::Array{Int64,2})

    avg_dps_all = Array{Float64}(undef,0)

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

    avg_dps_all = Array{Float64}(undef,0)

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

    low_dp_index_list = Array{Int64}(undef,0)

        for item = 1:length(sample_avg_list)
            if sample_avg_list[item] < 15
                push!(low_dp_index_list,item)
            end
        end

    low_dp_sample_names = Array{String}(undef,0)

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

    low_dp_index_list = Array{Int64}(undef,0)

        for item = 1:length(variant_avg_list)
            if variant_avg_list[item] < 15
                push!(low_dp_index_list,item)
            end
        end

    low_dp_positions = Array{Tuple{Any,Int64}}(undef,0)

        for i in low_dp_index_list
            chrom_position = chrom_labels[i,1],chrom_labels[i,2]
            push!(low_dp_positions,chrom_position)
        end

        return low_dp_positions



end

#functions for producing objects for plot functions

"""
    read_depth_threshhold(dp_array::Array{Int64,2})
Caps read depth outlier values at user defined threshhold.
threshhold defaults to dp = 100. All dp over 100 are set to 100 to visualize
read depth values between 0 < dp > 100 in better definition.
"""
function read_depth_threshhold(dp_array::Array{Int64,2})

    dp_array[dp_array[:,:].>100].=100

    return dp_array

end

"""
    save_numerical_array(num_array::Matrix{Any},sample_names,chr_labels,title,output_directory)
save numerical array with chr labels and sample ids to working directory
"""
function save_numerical_array(num_array,sample_names,chr_labels,title,output_directory)

      headings = hcat("chr,position")
      sample_names = hcat(headings,sample_names)

      chr_labeled_array_for_plotly=hcat(chr_labels, num_array)

      labeled_value_matrix_withsamplenames= vcat(sample_names,
                                            chr_labeled_array_for_plotly)

      writedlm(joinpath(output_directory ,"$(title)_matrix.csv"), labeled_value_matrix_withsamplenames, ",")
end

"""
    chromosome_label_generator(chromosome_labels::Array{Any,1})
Returns vector of chr labels and indices to mark chromosomes in plotly heatmap
Specifically, saves indexes and chrom labels in vectors to pass into heatmap
function to ticvals and tictext respectively.
Input is either gt_chromosome_labels or dp_chromosome_labels from
translate_gt/dp_to_num_array()
"""
function chromosome_label_generator(chromosome_labels::Array{Any,1})
    chrom_label_indices = findfirst.(map(a -> (y -> isequal(a, y)),
    unique(chromosome_labels)), [chromosome_labels])
    chrom_labels = unique(chromosome_labels)
    chrom_labels = [string(i) for i in chrom_labels]

    if length(chrom_labels) > 1
        for item=2:(length(chrom_labels))

            ratio=((chrom_label_indices[item])-(chrom_label_indices[item-1]))/(length(chromosome_labels))

            if ratio < 0.2
                font_size = "8"
                return chrom_labels,chrom_label_indices,font_size

            else
                font_size = "10"
                return chrom_labels,chrom_label_indices,font_size

            end
        end

    else

        font_size = "10"

        return chrom_labels,chrom_label_indices,font_size
    end
end

#=
"""
    chromosome_label_generator_for_chrom_bar(chromosome_labels::Array{String,2},number_rows)
Returns vector of chr labels and indices to build chromosome bar in heatmaply r
Specifically, saves indexes and chrom labels in vectors to pass into heatmap function to build array for row side colorbar
Input is either gt_chromosome_labels or dp_chromosome_labels from translate_gt/dp_to_num_array()
"""
function chromosome_label_generator_for_chrom_bar(chromosome_labels::Array{Any,1},number_rows)
    chrom_label_indices = findfirst.(map(a -> (y -> isequal(a, y)), unique(chromosome_labels)), [chromosome_labels])
    chrom_labels = unique(chromosome_labels)
    chrom_labels = [string(i) for i in chrom_labels]

    if chrom_label_indices[length(chrom_label_indices)] == number_rows
        duplicate_last_label="true"
    else
        duplicate_last_label="false"
    end

    push!(chrom_label_indices,number_rows)

    if length(chrom_labels) > 1
        for item=2:(length(chrom_labels))

            ratio=((chrom_label_indices[item])-(chrom_label_indices[item-1]))/(length(chromosome_labels))

            if ratio < 0.2
                font_size = "8"
                #println("font size is $font_size")
                return chrom_labels,chrom_label_indices,font_size,duplicate_last_label

            else
                font_size = "10"
                #println("font size is $font_size")
                return chrom_labels,chrom_label_indices,font_size,duplicate_last_label

            end
        end

    else

        font_size = "10"

        return chrom_labels,chrom_label_indices,font_size,duplicate_last_label
    end

end

=#

"""
    checkfor_outputdirectory(path::String)
Checks to see if output directory exists already. If it doesn't, it creates the
new directory to write output files to.
"""
function checkfor_outputdirectory(path::String)
    if isdir(path) == true
    else
            mkdir(path)
    end
end

"""
    generate_chromosome_positions_for_hover_labels(chr_labels::Array{Any,2})
creates tuple of genomic locations to set as tick labels. This is automatically
store chromosome positions in hover labels. However tick labels are set to hidden with showticklabels=false so they will not crowd the y axis.
"""
function generate_chromosome_positions_for_hover_labels(chr_labels::Array{Any,2})

returnXY_column1!(chr_labels)

chr_pos_tuple_list=Array{String}(undef,0)

    for row = 1:size(chr_labels,1)

        chr=chr_labels[row,1]
        pos=chr_labels[row,2]
        chr_pos_tuple="chr$chr,$pos"
        push!(chr_pos_tuple_list,chr_pos_tuple)
    end

    return chr_pos_tuple_list
end

#=function generate_chromosome_positions_for_hover_labels(chr_labels::Array{Any,2})

returnXY_column1!(chr_labels)

chr_pos_tuple_list=Array{Tuple}(undef,0)

    for row = 1:size(chr_labels,1)

        chr=chr_labels[row,1]
        pos=chr_labels[row,2]
        println(pos)
        chr_pos_tuple=chr,pos
        push!(chr_pos_tuple_list,chr_pos_tuple)
        println(chr_pos_tuple)
    end

    return chr_pos_tuple_list
end
=#

"""
    index_vcf(vcf_filename)
Creates and saves index file with three column array of vcf chrom, position, and row number to be used by significant list filter functions.
"""
function index_vcf(vcf_filename)

 reader = VCF.Reader(open(vcf_filename, "r"))

 chr_pos_col = Array{Any}(undef,0)
 row_col = Array{Any}(undef,0)

 row=1

  for record in reader

   chr = VCF.chrom(record)
   pos = VCF.pos(record)
   chr_pos_tuple = chr,pos
   row = row+1

   push!(chr_pos_col,chr_pos_tuple)
   push!(row_col,row)

  end

  index_matrix=hcat(chr_pos_col,row_col)
  return index_matrix

end

"""
    match_siglist_to_index(sig_list,vcf_index)
Returns vcf row indices of each variant position in sig_list for reader function to allow fast filtering in significant list filter funcitons.
"""
function match_siglist_to_index(sig_list,index_matrix)

 sig_list_tuples = Array{Any}(undef,0)

  for row = 1:size(sig_list,1)

   chrom = sig_list[row,1]
   pos = sig_list[row,2]
   sig_list_tuple = chrom,pos
   push!(sig_list_tuples,sig_list_tuple)

  end

  sig_list_index = findin(index_matrix,sig_list_tuples)

  return sig_list_index

end

"""
    make_chromosome_labels(chrom_label_info)
Returns vector of values to use as tick vals to show first chromosome label per chromosome with blank spaces between each first chromosome position for use with --y_axis_labels=chromosomes. duplicate_last_label tells if last chrom label is single or mutiple which affects number_to_fill value.
"""
function make_chromosome_labels(chrom_label_info)

       labels=chrom_label_info[1]
       index=chrom_label_info[2]
       duplicate_last_label=chrom_label_info[4]

       x=Array{Any}(undef,0)

       for n = 1:size(labels,1)

              push!(x, labels[n])
              counter=0
              number_to_fill=index[n+1]-index[n]

              while counter < number_to_fill - 1
                  counter = counter+1
                  push!(x, labels[n])
              end
        end

    if duplicate_last_label != "true"
      push!(x,labels[size(labels,1)])
    end

    return x

end

"""
    add_pheno_matrix_to_gt_data_for_plotting(pheno_matrix,gt_num_array,trait_labels,chrom_label_info,number_rows)
add the pheno matrix used to group samples to the data array for input into plotting functions. Resizes the pheno matrix to maintain correct dimensions for heatmap viz by finding value=0.05*number_rows_data to multiply each pheno row by before vcat.
"""
function add_pheno_matrix_to_gt_data_for_plotting(gt_num_array,pheno_matrix,trait_labels,chrom_label_info,number_rows)

    pheno_matrix=pheno_matrix[2:size(pheno_matrix,1),:]

    #consider options when have few variants. no multiplyer when can't find even multiple?

    #find multiplyer value
    pheno_row_multiplyer=0.05*(size(gt_num_array,1))
    pheno_row_multiplyer=(pheno_row_multiplyer)
    pheno_row_multiplyer=round(pheno_row_multiplyer/(size(pheno_matrix,1)))
    #println(size(pheno_matrix,1))

    for i = 1:(size(pheno_matrix, 1))
        for n=1:(size(pheno_matrix, 2))
            pheno_matrix[i, n] = pheno_matrix[i, n] == 1 ? -1 : pheno_matrix[i, n]
            pheno_matrix[i, n] = pheno_matrix[i, n] == 2 ? -2 : pheno_matrix[i, n]
        end
    end

    #resize pheno matrix and create trait_label_array for labeling
    resized_pheno_matrix=pheno_matrix[1:1, :]

    trait_label_array=trait_labels[1]

    for row = 1:size(pheno_matrix,1)

        trait_row=pheno_matrix[row:row, :]
        trait = trait_labels[row]

        i=0

        while i <= (pheno_row_multiplyer+1)
            i=i+1
            #println(i)
            resized_pheno_matrix=vcat(resized_pheno_matrix,trait_row)
            trait_label_array=vcat(trait_label_array,trait)
        end
    end

    trait_label_indices = findfirst.(map(a -> (y -> isequal(a, y)),
    unique(trait_label_array)), [trait_label_array])
    pheno_trait_labels = unique(trait_label_array)
    pheno_trait_labels = [String(i) for i in pheno_trait_labels]
    trait_label_indices = [i+number_rows for i in trait_label_indices]

    chrom_labels=chrom_label_info[1]
    chrom_label_indices=chrom_label_info[2]

    chrom_labels = vcat(chrom_labels,pheno_trait_labels)
    chrom_label_indices = vcat(chrom_label_indices,trait_label_indices)

    chrom_label_info = chrom_labels,chrom_label_indices,chrom_label_info[3]

    data_and_pheno_matrix=vcat(gt_num_array,resized_pheno_matrix)

    data_and_pheno_matrix=convert(Array{Int64,2},data_and_pheno_matrix)

    return data_and_pheno_matrix,trait_label_array,chrom_label_info
end

"""
    add_pheno_matrix_to_dp_data_for_plotting(pheno_matrix,dp_num_array,trait_labels,chrom_label_info,number_rows)
add the pheno matrix used to group samples to the data array for input into plotting functions. Resizes the pheno matrix to maintain correct dimensions for heatmap viz by finding value=0.05*number_rows_data to multiply each pheno row by before vcat.
"""
function add_pheno_matrix_to_dp_data_for_plotting(dp_num_array,pheno_matrix,trait_labels,chrom_label_info,number_rows)

    pheno_matrix=pheno_matrix[2:size(pheno_matrix,1),:]

    #consider options when have few variants. no multiplyer when can't find even multiple?

    #find multiplyer value
    pheno_row_multiplyer=0.05*(size(dp_num_array,1))
    pheno_row_multiplyer=round(pheno_row_multiplyer)
    pheno_row_multiplyer=pheno_row_multiplyer/(size(pheno_matrix,1))

    for i = 1:(size(pheno_matrix, 1))
        for n=1:(size(pheno_matrix, 2))
            pheno_matrix[i, n] = pheno_matrix[i, n] == 1 ? -40 : pheno_matrix[i, n]
            pheno_matrix[i, n] = pheno_matrix[i, n] == 2 ? -60 : pheno_matrix[i, n]
        end
    end

    #resize pheno matrix and create trait_label_array for labeling
    resized_pheno_matrix=pheno_matrix[1:1, :]

    trait_label_array=trait_labels[1]

    for row = 1:size(pheno_matrix,1)

        trait_row=pheno_matrix[row:row, :]
        trait = trait_labels[row]

        i=0

        while i <= (pheno_row_multiplyer+1)
            i=i+1
            resized_pheno_matrix=vcat(resized_pheno_matrix,trait_row)
            trait_label_array=vcat(trait_label_array,trait)
        end
    end

    trait_label_indices = findfirst.(map(a -> (y -> isequal(a, y)),
    unique(trait_label_array)), [trait_label_array])
    pheno_trait_labels = unique(trait_label_array)
    pheno_trait_labels = [String(i) for i in pheno_trait_labels]
    trait_label_indices = [i+number_rows for i in trait_label_indices]

    chrom_labels=chrom_label_info[1]
    chrom_label_indices=chrom_label_info[2]

    chrom_labels = vcat(chrom_labels,pheno_trait_labels)
    chrom_label_indices = vcat(chrom_label_indices,trait_label_indices)

    chrom_label_info = chrom_labels,chrom_label_indices,chrom_label_info[3]

    #vcat pheno matrix to gt_num_array
    data_and_pheno_matrix=vcat(dp_num_array,resized_pheno_matrix)

    data_and_pheno_matrix=convert(Array{Int64,2},data_and_pheno_matrix)

    return data_and_pheno_matrix,trait_label_array,chrom_label_info
end

"""
    generate_hover_text_array(chr_pos_tuple_list,sample_names,input,mode)
Generate array of data for hovertext to use as custom hover text for ungrouped heatmaps. Where mode is GT or DP.
"""
function generate_hover_text_array(chr_pos_tuple_list,sample_names,input,mode)

hover_text_array=Array{Any,1}(undef,0)

    for n=1:length(chr_pos_tuple_list)
        pos=chr_pos_tuple_list[n]

        row_data_vector=Array{Any,1}(undef,0)

        for m=1:length(sample_names)
            id=sample_names[m]

            z=input[n,m]

            if mode == "GT"
            z_dict = Dict(0 => "No Call", 1 => "Homozygous Reference", 2 => "Heterozygous Variant", 3 => "Homozygous Variant")
            z=z_dict[z]

            elseif mode == "DP"
                if z==-1
                    z="No Call"
                end
            end

            text_cell="ID: $id \n $pos \n $mode: $z"

            push!(row_data_vector,text_cell)
        end

        push!(hover_text_array,row_data_vector)
    end

return(hover_text_array)

end

"""
    generate_hover_text_array_grouped(chr_pos_tuple_list,sample_names,input,mode)
Generate array of data for hovertext to use as custom hover text for grouped heatmaps. Where mode is GT or DP.
"""
function generate_hover_text_array_grouped(chr_pos_tuple_list,sample_names,input,mode,number_rows)

hover_text_array=Array{Any,1}(undef,0)

    for n=1:number_rows
        pos=chr_pos_tuple_list[n]

        row_data_vector=Array{Any,1}(undef,0)

        for m=1:length(sample_names)
            id=sample_names[m]

            z=input[n,m]

            if mode == "GT"
            z_dict = Dict(0 => "No Call", 1 => "Homozygous Reference", 2 => "Heterozygous Variant", 3 => "Homozygous Variant")
            z=z_dict[z]

            elseif mode == "DP"
                if z==-1
                    z="No Call"
                end
            end

            text_cell="ID: $id \n $pos \n $mode: $z"

            push!(row_data_vector,text_cell)
        end

        push!(hover_text_array,row_data_vector)
    end

    for n=number_rows+1:length(chr_pos_tuple_list)
        trait=chr_pos_tuple_list[n]

        split_trait=split(trait,',')

        row_data_vector=Array{Any,1}(undef,0)

        for m=1:length(sample_names)
            id=sample_names[m]

            z=input[n,m]

            if mode == "GT"
                if z==-1
                    z="Trait 1: $(split_trait[1])"
                elseif z==-2
                    z="Trait 2: $(split_trait[2])"
                end

            elseif mode == "DP"
                if z==-40
                    z="Trait 1: $(split_trait[1])"
                elseif z==-60
                    z="Trait 2: $(split_trait[2])"
                end
            end

            text_cell="ID: $id \n $z"

            push!(row_data_vector,text_cell)
        end

        push!(hover_text_array,row_data_vector)
    end

return(hover_text_array)

end

"""
    save_graphic(graphic,output_directory,save_ext,title,remote_option)
Save plot in either html or static image formats incuding eps, png, svg, and pdf
"""
function save_graphic(graphic,output_directory,save_ext,title,remote_option)

    if save_ext=="html" && remote_option == true
        PlotlyJS.savehtml(graphic, joinpath("$(output_directory)" ,"$title.$(save_ext)"), :remote)
    elseif save_ext=="html" && remote_option != true
        PlotlyJS.savehtml(graphic, joinpath("$(output_directory)" ,"$title.$(save_ext)"))
    elseif save_ext != "html"
        PlotlyJS.savefig(graphic, joinpath("$(output_directory)" ,"$title.$(save_ext)"))
    end

    #display(graphic)

end

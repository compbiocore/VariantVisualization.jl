#guide for writing technical docs - https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php

# *** Outline of program: ***

#A.) load and clean VCF

             #1) determine row number of header line to avoid "##" in META intro lines, read vcf file, clean and order in ascending chr position order
             #2) data cleaning

#B) define all functions

        #1) field selection

              #a) FORMAT reader - get index of fields to visualize
              #b) maf match list for maf correction of genotype
              #c) genotype selection with maf correction of genotype
              #d) genotype selection with no maf correction
              #e) read depth selection

        #2) variant selection

              #a) siglist match filter
              #b) match chromosome range filter
              #c) rare variants only

        #3) rearrange and select columns

              #a) rearrange columns by phenotype with use input phenotype matrix
              #b) select columns to visualize

        #4) plotting functions

              #a) define plotlyJS function for genotype heatmap
              #b) define plotlyJS function for read depth heatmap

#C) run functions for every possible combination of features

        #1) universal filters

               #a) PASS FILTER variants only
               #b) reorder columns to match list
               #c) select columns to visualize

        #2) combinations of field selection and variant selection filters

               #a) genotype / display all variants
               #b) genotype / select variants matching significant variant location list
               #c) genotype / select variants within chromosomal range
               #d) read depth / display all variants
               #e) read depth / select variants matching significant variant location list
               #f) read depth / select variants within chromosomal range

# *** End of outline ***

using DataFrames #use CSV.jl ? depwarnings
using PlotlyJS
using Rsvg
using Blink

g_white = "400" #homo reference 0/0
g_red = "800" #homo variant 1/1 1/2 2/2 1/3 2/3 3/3 4/4 5/5 6/6 etc
g_pink = "600" #hetero variant 0/1 1/0 0/2 2/0 etc
g_blue = "0" #no call ./.

#A.) load and clean VCF

#1) determine row number of header line to avoid "##" in META intro lines, read vcf file, clean and order in ascending chr position order

readvcf = readlines(ARGS[1])

    for row = 1:size(readvcf,1)

        if contains(readvcf[row], "#CHROM")
            header_col = row
            global header_col
            header_string = readvcf[row]
        end
    end

skipstart_number=header_col-1 #this allows vcf to be loaded by readtable, the META lines at beginning of file start with "##" and interfere with readtable function - need readtable vs readdlm for reorder columns
df_vcf=readtable(ARGS[1], skipstart=skipstart_number, separator='\t')

vcf=Matrix(df_vcf)

#load vcf as dataframe twice: once to use for matrix and other to pull header info from

#df_vcf=readtable(ARGS[1], skipstart=skipstart_number, separator='\t')

#2) data cleaning

    #replace String "X" and "Y" from chromosome column so all elements are Int64 type
    function clean_column1!(x)
        n = size(x, 1)
        for i = 1:n
            x[i, 1] = x[i, 1] == "X" ? 23 : x[i, 1]
            x[i, 1] = x[i, 1] == "Y" ? 24 : x[i, 1]
        end
    end

clean_column1!(vcf) #cleancolumn!() replaces "X" and "Y" with Int64 23 and 24, rest of chr numbers are not loaded as Int64 and need to be same type

    for n = 1:size(vcf,1)
        if vcf[n, 1] != 23
        vcf[n, 1] = parse(Int64, vcf[n, 1])
        end
    end
#=
    for n = 1:size(vcf,1)
        if vcf[n, 1] != 24
        vcf[n, 1] = parse(Int64,vcf[n, 1])
        end
    end
=#
    #sort rows by chr then chromosome position so are in order of chromosomal architecture
vcf = sortrows(vcf, by=x->(x[1],x[2]))

#B) define all functions

#1) field selection

    #a) FORMAT reader - get index of fields to visualize

    #what features to visualize from INFO field - or others - if FORMAT / genotype info isn't included
    #index other data per variant for bars
    #check if format col is included, if so - run function, if not are there cells for each sample? I don't think so - check this.

    function format_reader()
        format = vcf[1,9]

        s = split(format,":")
        gt_index = find(x -> x == "GT",s)
        dp_index = find(x -> x == "DP",s)
        gq_index = find(x -> x == "GQ",s) #genotype quality - look at annotated pdf to determine how to interpret
        pl_index = find(x -> x == "PL",s)
        mq_index = find(x -> x == "MQ",s)

            if ARGS[3] == "-gt"
                index = gt_index
            elseif ARGS[3] == "-dp"
                index = dp_index
            elseif ARGS[3] == "-gq"
                index = gq_index
            elseif ARGS[3] == "-pl"
                index = pl_index
            elseif ARGS[3] == "-mq"
                index = mq_index
            end

        return index
    end

    index = format_reader()

        #ad not necessary :
        #=DP is equal to the number of reads used for calling the variant,
        AD numbers are the total number of reads with each allele (including reads
        not used for calling due to low quality or whatever).=#

    #b) maf match list for maf correction of genotype

#=
    vcf = DataFrame(vcf)
    println("Loading maf list...")
    #maf_list=readtable("combined_EXaC_1000G_maf_list.txt",separator = ' ') - takes LONG time
    maf_list=readdlm("combined_EXaC_1000G_maf_list.txt")
    maf_list = DataFrame(maf_list)
    println("Finished loading maf list.")

    @time match_list = join(maf_list, vcf, on = [:x1, :x2], kind = :right)
=#

#c) genotype selection with maf correction of genotype

    #may need to define vcf first, or get rid of conditional eval and set 'maf_sub = maf_list_match_vcf(vcf)' inside function definition

    function genotype_cell_searcher_maf_correction(x) #when x is output subarray of variant selection

        for row = 1:size(x,1)

            maf = match_list[row,3]

                for col = 10:size(x,2)

                    cell_to_search = x[row,col]
                    S = split(cell_to_search, ":")
                    genotype = S[index]

                    homo_variant = ["1/1" "1/2" "2/2" "1/3" "2/3" "3/3" "1/4" "2/4" "3/4" "4/4" "1/5" "2/5" "3/5" "4/5" "5/5" "1/6" "2/6" "3/6" "4/6" "5/6" "6/6" "1|1" "1|2" "2|2" "1|3" "2|3" "3|3" "1|4" "2|4" "3|4" "4|4" "1|5" "2|5" "3|5" "4|5" "5|5" "1|6" "2|6" "3|6" "4|6" "5|6" "6|6"]
                    hetero_variant = ["0/1" "0/2" "0/3" "0/4" "0/5" "0/6" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "0|1" "0|2" "0|3" "0|4" "0|5" "0|6" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0"]

                    if genotype == "./." || ".|."
                        x[row,col] = g_blue

                    elseif genotype == "0/0" || "0|0"
                        x[row,col] = g_white

                    else
                        for i in hetero_variant
                            if genotype == i
                                x[row,col] = g_pink
                            end
                        end

                        for i in homo_variant
                            if genotype == i
                                if maf > (0.5) && maf < (1.1)|| MAF==(0.5)
                                    x[row,col] = g_white
                                elseif maf < 0.5 && MAF >= 0
                                    x[row,col] = g_red
                                end
                            end
                        end

                    end
                end
        end
        return x
    end

    #d) genotype selection with no maf correction

    function genotype_cell_searcher(x) #when x is output subarray of variant selection

        for row=1:size(x,1)

            for col = 10:size(x,2)

                cell_to_search = x[row,col]
                S = split(cell_to_search, ":")
                genotype = S[index[1]]

                homo_variant = ["1/1" "1/2" "2/2" "1/3" "2/3" "3/3" "1/4" "2/4" "3/4" "4/4" "1/5" "2/5" "3/5" "4/5" "5/5" "1/6" "2/6" "3/6" "4/6" "5/6" "6/6" "1|1" "1|2" "2|2" "1|3" "2|3" "3|3" "1|4" "2|4" "3|4" "4|4" "1|5" "2|5" "3|5" "4|5" "5|5" "1|6" "2|6" "3|6" "4|6" "5|6" "6|6"]
                hetero_variant = ["0/1" "0/2" "0/3" "0/4" "0/5" "0/6" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "0|1" "0|2" "0|3" "0|4" "0|5" "0|6" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0"]

                if genotype == "./." #|| ".|."
                    x[row,col] = g_blue

                elseif genotype == "0/0" #|| "0|0"
                    x[row,col] = g_white

                else

                    for i in hetero_variant
                        if genotype == i
                            x[row,col] = g_pink
                        end
                    end

                    for i in homo_variant
                        if genotype == i
                            x[row,col] = g_red
                        end
                    end
                end
            end
        end
        return x
    end

    #e) read depth selection maf correction

    function dp_cell_searcher(x) #when x is output subarray of variant selection

        for row=1:size(x,1)

            for col=10:size(x,2)

                dp_cell=x[row,col]
                S=split(dp_cell, ":")
                dp=S[3]

                    if x[row,col] == x[row,col]
                        x[row,col]=dp
                    end

            end
        end

        return x
    end

#2) variant selection

    #a) siglist match filter
    function sig_list_vcf_filter(y,x) #where x = significant list ordered with 23s y is vcf - y is vcf

        sig_list_subarray_pre=Array{Any}(0,(size(y,2)))

            for row= 1:size(x,1)

                chr=x[row,1]
                pos=x[row,2]
                first_sub=y[(y[:,1].== chr),:]
                final_sub=first_sub[first_sub[:,2].== pos,:]
                sig_list_subarray_pre = [sig_list_subarray_pre; final_sub]
            end

            return sig_list_subarray_pre
    end

#= oringial siglist match - need to fix to mae modular
#a) siglist match filter
function sig_list_vcf_filter(y,x) #where x = significant list ordered with 23s y is vcf - y is vcf

    sig_list_subarray_pre=Array{Any}(0,(size(vcf,2)))

        df1_vcf=DataFrame(vcf)
        println(df1_vcf)
        for row= 1:size(x,1)

            chr=x[row,1]
            pos=x[row,2]
            first_sub=vcf[(vcf[:,1].== chr),:]
            final_sub=first_sub[first_sub[:,2].== pos,:]
            sig_list_subarray_pre = [sig_list_subarray_pre; final_sub]
        end

        return sig_list_subarray_pre
end

=#

    #b) match chromosome range filter
    function chromosome_range_vcf_filter(x::AbstractString)

            a=split(x,":")
            chrwhole=a[1]
            chrnumber=split(chrwhole,"r")
            string_chr=chrnumber[2]
            chr=parse(string_chr)
            range=a[2]
            splitrange=split(range, "-")
            lower_limit=splitrange[1]
            lower_limit_int=parse(lower_limit)
            upper_limit=splitrange[2]
            upperlimit_int=parse(upper_limit)

            #make subarray within range

            first_sub=vcf[(vcf[:,1].== chr),:]
            second_sub=first_sub[(first_sub[:,2].>= lower_limit_int),:]
            chr_range_subarray=second_sub[(second_sub[:,2].<= upperlimit_int),:]
            return chr_range_subarray

    end

    #c) rare variants only

    #=
        function rare_variant_vcf_filter(x)

            rare_variant_subarray_pre=Array{Any}(0,(size(vcf,2)))

            for row = 1:size(sub_array_maf)
                maf = sub_array_maf[row,3]
                rare_variant_subarray = vcf[(sub_array_maf[:,3].< 0.01),:]
            end
    =#

#3) rearrange and select columns

    #if add feature to select columns, do this before filtering rows to make filter functions perform faster? (fewer columns)

#a) rearrange columns by phenotype with use input phenotype matrix

    #phenotype matrix - format must be CSV, 0's will be sorted first

    function load_sort_phenotype_matrix(x,y) #when x is ARGS[6] which is phenotype_matrix which is in *CSV format! and y is pheno row to sort on

        global vcf
        global df_vcf

        pheno = readdlm(x, ',')

        #get row numer of user-chosen phenotype characteristic to sort columns by
        row_to_sort_by = find(x -> x == y, pheno)
        row_to_sort_by = row_to_sort_by[1]

        #remove phenotype_row_labels used to identify row to sort by, so row can be sorted without strings causing issues
        pheno = pheno[:,2:size(pheno,2)]

        pheno = sortcols(pheno, by = x -> x[row_to_sort_by], rev = false)

        id_list = pheno[1,:]
        vcf_header = names(df_vcf)
        vcf_info_columns = vcf_header[1:9]

        for item = 1:size(id_list,2) #names in sample list dont match vcf so have to clean
               id_list[item] = replace(id_list[item],".","_")
               id_list[item] = Symbol(id_list[item])
        end

        new_order = vcat(vcf_info_columns,id_list)

        header_as_strings = new_order

        col_new_order=vec(new_order)

        df1_vcf = DataFrame(vcf)

        rename!(df1_vcf, f => t for (f, t) = zip(names(df1_vcf), names(df_vcf)))

        vcf = df1_vcf[:, col_new_order]

        vcf = Matrix(vcf)

        return vcf
    end


    function reorder_columns(x) #where x = ARGS[#] which is name of new order list
        global vcf
        global df_vcf

        new_order=readdlm(x)

        header_as_strings = new_order

        for item = 1:size(new_order,2) #names in sample list dont match vcf so have to clean
               new_order[item] = replace(new_order[item],".","_")
               new_order[item] = Symbol(new_order[item])
        end

        col_new_order=vec(new_order)

        df1_vcf = DataFrame(vcf)

        rename!(df1_vcf, f => t for (f, t) = zip(names(df1_vcf), names(df_vcf)))

        vcf = df1_vcf[:, col_new_order]

        vcf = Matrix(vcf)

        return vcf
    end

#b) select columns to visualize

function select_columns(x) #where x = ARGS[#] which is name of list of samples to include

    global vcf
    global df_vcf

    selectedcolumns=readdlm(x)

    header_as_strings = selectedcolumns

    for item = 1:size(selectedcolumns,2) #names in sample list dont match vcf so have to clean
           selectedcolumns[item] = replace(selectedcolumns[item],".","_")
           selectedcolumns[item] = Symbol(selectedcolumns[item])
    end

    col_selectedcolumns=vec(selectedcolumns)

    df1_vcf = DataFrame(vcf)

    rename!(df1_vcf, f => t for (f, t) = zip(names(df1_vcf), names(df_vcf)))

    vcf = df1_vcf[:, col_selectedcolumns]

    vcf = Matrix(vcf)

    return vcf
end


#4) plotting functions

    #a) define plotlyJS function for genotype heatmap

    function genotype_heatmap2(x) #when x = array_for_plotly

        trace=heatmap(
            z = x,
            transpose=true,
            #colorscale = "Picnic",
            colorscale = "Rainbow",
            gridcolor = "#E2E2E2",
            showscale = true
            );

        layout = Layout(
                        title = "$title",#defines title of plot
                        xaxis=attr(title="Sample Number", showgrid=false, zeroline=false),
                        yaxis=attr(title="Chromosomal Location", zeroline=false)
        )

        data = (trace)
        plot(data,layout) #call plot type and layout with all attributes to plot function
    end

    #b) define plotlyJS function for read depth heatmap

    function dp_heatmap2(x) #when x = array_for_plotly

        trace=heatmap(
            z = x,
            transpose=true,
            colorscale = "YIGnBu",
            gridcolor = "#E2E2E2",
            showscale = true
            );

        layout = Layout(
                        title = "$title",#defines title of plot
                        xaxis=attr(title="Sample Number", showgrid=false, zeroline=false),
                        yaxis=attr(title="Chromosomal Location", zeroline=false)
        )

        data = (trace)
        plot(data,layout) #call plot type and layout with all attributes to plot function
    end

#C) Run functions for every possible combination of features

#1) universal filters

    #a) PASS FILTER variants only

    if ARGS[5] == "pass_only"
        vcf=vcf[(vcf[:,7].== "PASS"),:]
    end

    #b) reorder columns to match list

    if ARGS[6] == "reorder_columns"
        #vcf = reorder_columns(ARGS[7])
        vcf = load_sort_phenotype_matrix(ARGS[7],ARGS[11])
    end

    #c) select columns to visualize

    if ARGS[9] == "select_columns"
        vcf = select_columns(ARGS[10])
    end


#2) combinations of field selection and variant selection filters

    #a) genotype / display all variants

    if ARGS[3] == "-gt" && ARGS[4] == "-a"

        #replace cells of vcf file with representative values for field chosen (genotype value)
        vcf = genotype_cell_searcher(vcf)

        #convert value overwritten vcf into subarray of just values, no annotation/meta info
        array_for_plotly=vcf[:,10:size(vcf,2)]

        #define title for plot
        title = "Genotype Data for All Variants"

        #plot heatmap for genotype and save as format specified by ARGS[2], defaults to pdf
        graphic = genotype_heatmap2(array_for_plotly)
        extension=ARGS[2] #must define this variable, if use ARGS[2] directly in savefig it is read as String[pdf] or something instead of just "pdf"
        PlotlyJS.savefig(graphic, "all_genotype.$extension")

        #=***
        activate this block if want to export labeled value matrix for internal team use

        df_withsamplenames=readtable(ARGS[1], skipstart=skipstart_number, header=false,separator='\t')
        samplenames=df_withsamplenames[1,10:size(df_withsamplenames,2)]
        samplenames=Matrix(samplenames)
        chrlabels=vcf[:,1:2]

        chr_labeled_array_for_plotly=hcat(chrlabels, array_for_plotly)
        labeled_value_matrix_withsamplenames= vcat(samplenames,chr_labeled_array_for_plotly)

        writedlm("labeled_value_matrix.txt", labeled_value_matrix_withsamplenames, "\t")
        ***=#

    elseif ARGS[3] == "-gt" && ARGS[4] == "-l"
        df1=DataFrame(vcf)

        #load siglist file
        siglist_unsorted=readdlm(ARGS[8], ',',skipstart=1)

        #replace X with 23 and sort by chr# and position
        clean_column1!(siglist_unsorted)
        siglist=sortrows(siglist_unsorted, by=x->(x[1],x[2]))

        #significant variants only filter - subarray of vcf matching lsit of variants of interest
        vcf = sig_list_vcf_filter(vcf,siglist)

        #write over vcf to create value matrix for genotype fieldtype
        sig_list_subarray_post=genotype_cell_searcher(vcf)

        #convert value overwritten vcf into subarray of just values, no annotation/meta info
        array_for_plotly=sig_list_subarray_post[:,10:size(sig_list_subarray_post,2)]
        title = "Genotype Data for Variants of Interest"

        #plot heatmap for genotype and save as format specified by ARGS[2], defaults to pdf
        graphic = genotype_heatmap2(array_for_plotly)
        extension=ARGS[2] #must define this variable, if use ARGS[5] directly in savefig it is read as String[pdf] or something instead of just "pdf"
        PlotlyJS.savefig(graphic, "siglist_genotype.$extension")

    elseif ARGS[3] == "-gt" && ARGS[4] == "-r"

        #define range of variants to visualize
        chr_range = ARGS[8]

        #create subarray of vcf matching range parameters
        chr_range_subarray_pre = chromosome_range_vcf_filter(chr_range)

        #write over vcf to create keyed-values matrix showing genotype
        chr_range_subarray_post = genotype_cell_searcher(chr_range_subarray_pre)

        #convert value overwritten vcf into subarray of just values, no annotation/meta info
        array_for_plotly=chr_range_subarray_post[:,10:size(chr_range_subarray_post,2)]

        chrlabels=chr_range_subarray_post[:,1:2]
        chr_labeled_array_for_plotly=hcat(chrlabels, array_for_plotly)      #array for R_translation
        writedlm("labeled_value_matrix_chr_range.txt", chr_labeled_array_for_plotly, "\t")

        #define title
        title = "Genotype Data for Variants within $(ARGS[8])"

        #plot heatmap for genotype and save as format specified by ARGS[2], defaults to pdf
        graphic = genotype_heatmap2(array_for_plotly)
        extension=ARGS[2] #must define this variable, if use ARGS[2] directly in savefig it is read as String[pdf] or something instead of just "pdf"
        PlotlyJS.savefig(graphic, "chr_range_genotype.$extension")

    elseif ARGS[3] == "-dp" && ARGS[4] == "-a"

        #replace cells of vcf file with representative values for field chosen (genotype value)
        vcf = dp_cell_searcher(vcf)

        #convert value overwritten vcf into subarray of just values, no annotation/meta info
        array_for_plotly=vcf[:,10:size(vcf,2)]

        #define title for plot
        title = "Read Depth Data for All Variants"

        #plot heatmap for genotype and save as format specified by ARGS[2], defaults to pdf
        graphic = dp_heatmap2(array_for_plotly)
        extension=ARGS[2] #must define this variable, if use ARGS[2] directly in savefig it is read as String[pdf] or something instead of just "pdf"
        PlotlyJS.savefig(graphic, "all_readdepth.$extension")

        #=***
        activate this block if want to export labeled value matrix for internal team use

        df_withsamplenames=readtable(ARGS[1], skipstart=skipstart_number, header=false,separator='\t')
        samplenames=df_withsamplenames[1,10:size(df_withsamplenames,2)]
        samplenames=Matrix(samplenames)
        chrlabels=vcf[:,1:2]

        chr_labeled_array_for_plotly=hcat(chrlabels, array_for_plotly)
        labeled_value_matrix_withsamplenames= vcat(samplenames,chr_labeled_array_for_plotly)

        writedlm("labeled_value_matrix.txt", labeled_value_matrix_withsamplenames, "\t")
        ***=#

    elseif ARGS[3] == "-dp" && ARGS[4] == "-l"

        #load siglist file
        siglist_unsorted=readdlm(ARGS[8], ',',skipstart=1)

        #replace X with 23 and sort by chr# and position
        clean_column1!(siglist_unsorted)
        siglist=sortrows(siglist_unsorted, by=x->(x[1],x[2]))

        #create subarray of vcf per siglist
        sig_list_subarray_pre=sig_list_vcf_filter(siglist)

        #write over vcf to create keyed-values matrix showing genotype
        sig_list_subarray_post=dp_cell_searcher(sig_list_subarray_pre)

        #convert value overwritten vcf into subarray of just values, no annotation/meta info
        array_for_plotly=sig_list_subarray_post[:,10:size(sig_list_subarray_post,2)]
        title = "Read Depth Data for Variants of Interest"

        #plot heatmap for read depth and save as format specified by ARGS[2], defaults to pdf
        graphic = dp_heatmap2(array_for_plotly)
        extension=ARGS[2] #must define this variable, if use ARGS[2] directly in savefig it is read as String[pdf] or something instead of just "pdf"
        PlotlyJS.savefig(graphic, "siglist_readdepth.$extension")

    elseif ARGS[3] == "-dp" && ARGS[4] == "-r"

        #define range of variants to visualize
        chr_range = ARGS[8]

        #create subarray of vcf matching range parameters
        chr_range_subarray_pre = chromosome_range_vcf_filter(chr_range)

        #write over vcf to create keyed-values matrix showing genotype
        chr_range_subarray_post = dp_cell_searcher(chr_range_subarray_pre)

        #convert value overwritten vcf into subarray of just values, no annotation/meta info
        array_for_plotly=chr_range_subarray_post[:,10:size(chr_range_subarray_post,2)]

        #define title
        title = "Read Depth Data for Variants within $(ARGS[8])"

        #plot heatmap for read depth and save as format specified by ARGS[2], defaults to pdf
        graphic = dp_heatmap2(array_for_plotly)
        extension=ARGS[2] #must define this variable, if use ARGS[5] directly in savefig it is read as String[pdf] or something instead of just "pdf"
        PlotlyJS.savefig(graphic, "chr_range_readdepth.$extension")

        #Conditional Genotype quality (GQ) ex. 12

        #Phred-scaled genotype likelihood (PL) ex. 79

        #RMS mapping quality (MQ)


    end

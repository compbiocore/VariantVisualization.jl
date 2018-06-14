
function variant_filter!(variant_filter::String, vcf::Matrix, list="", chr_range="")
    if variant_filter == "pass_only"
        vcf=vcf[(vcf[:,7].== "PASS"),:]
    elseif variant_filter == "list"
        #function to filter vcf with list
    elseif variant_filter == "range"
        if chr_range == ""
            println("Enter a chromosome range")
        else
            chr_range = #regex to confirm ch1:20000000-30000000
            println("selecting variants within $chr_range")
            vcf = ViVa.chromosome_range_vcf_filter(chr_range,vcf)
        end

    end
    return vcf
end


"""
jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)

filters, plots visualization, and saves as figure.
utilizes all global variables set in first cell of jupyter notebook

"""

function jupyter_main(vcf_filename,field_to_visualize,variant_filter,sample_filter,save_format,plot_title)

#A)Load vcf and identify field index
vcf_tuple = ViVa.load_vcf(vcf_filename)
original_vcf = vcf_tuple[1]
df_vcf = vcf_tuple[2]

index = ViVa.format_reader(original_vcf, field_to_visualize)

#retain original version of vcf for second run
vcf = original_vcf

#B)Apply filters # issue - if vector need to use in, if string need to use contains() / sample_filters may be both - how to conditionally use correct function

#1)Variant filters

#=
for i = 1:size(variant_filter,1)
    vcf = filter(variant_filter[i], vcf)
end
=#

if typeof(variant_filter) != "String"

for i = 1:size(variant_filter,1)

    if variant_filter[i] == "pass_only"
        println("selecting pass_only variants")
        vcf=vcf[(vcf[:,7].== "PASS"),:]

    elseif variant_filter[i] == "range"
        chr_range = variant_filter[i+1]
        println("selecting variants within $chr_range")
        vcf = ViVa.chromosome_range_vcf_filter(chr_range,vcf)

    elseif variant_filter[i] == "list"
        siglist_file = variant_filter[i+1]

        function load_siglist(x)

        siglist_unsorted = readdlm(siglist_file, ',', skipstart=1)
        ViVa.clean_column1!(siglist_unsorted)
        siglist = sortrows(siglist_unsorted, by=x->(x[1],x[2]))

        return siglist
    end

    siglist = (load_siglist(siglist_file))

    vcf = ViVa.sig_list_vcf_filter(vcf,siglist)

end
end

else 
    if variant_filter == "pass_only"
        println("selecting pass_only variants")
        vcf=vcf[(vcf[:,7].== "PASS"),:]

    elseif variant_filter == "range"
        println("you must enter a chromosome range in format chr1:20000000-30000000")

    elseif variant_filter[i] == "list"
        println("you must enter a list filename in tab delimited format")

    end
end

#2) Sample Filters
for i = 1:size(sample_filter,1)
    if sample_filter[i] == "reorder_columns"
        list =sample_filter[i+1]
        key = sample_filter[i+2]
        println("sorting columns by $key in $list")
        vcf = load_sort_phenotype_matrix(list, key, vcf, df_vcf)

    elseif sample_filter[i] == "select_columns"
        id_list = sample_filter[i+1]
        println("selecting columns to match $id_list")
        vcf = ViVa.select_columns(id_list, vcf, df_vcf)
    end
end

#C) Convert filtered vcf to value matrix then plot and save
if field_to_visualize == "genotype"

    value_matrix = ViVa.genotype_cell_searcher(vcf,index)
    array_for_plotly=value_matrix[:,10:size(value_matrix,2)]
    graphic = ViVa.genotype_heatmap2(array_for_plotly,plot_title)
    PlotlyJS.savefig(graphic, "$plot_title.$save_format")

elseif field_to_visualize == "read_depth"

    value_matrix=ViVa.dp_cell_searcher(vcf,index)
    array_for_plotly=value_matrix[:,10:size(value_matrix,2)]
    graphic = ViVa.dp_heatmap2(array_for_plotly,plot_title)
    PlotlyJS.savefig(graphic, "$plot_title.$save_format")

end

end #end of jupyter_main()







#check to see if is a vector or string
function isvector(var)

 if typeof(var) == String
  println("String")

 else
  number_inputs = (length(var))
  println(number_inputs)
 end
end



#check that variables have all necessary arguments
#=
if contains("select_columns",sample_filter) # must have "sample_id_list.txt", "reorder_columns", must have "phenotype_matrix.csv" & "phenotype_label_to_sort_by
 if contains("list.txt", sample_filter)
 else
  println("you must enter a list of sample id's to select for visualization")
 end
end

elseif Tuple in (typeof(sample_filter))
 println("tuple")
else
 println("keep trying")
end


if "reorder_columns" in sample_filter
 if "phenotype_matrix.csv" in sample_filter
  if "phenotype_label_to_sort_by" in sample_filter
  else
   println("you must enter a phenotype_label_to_sort_by")
  end
 else
  println("phenotype_matrix.csv")
 end
end


 else println("you must enter a list of sample id's to select for visualization")
 end
end

=#

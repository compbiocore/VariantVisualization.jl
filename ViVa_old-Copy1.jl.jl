
vcf_filename = "variants.filtered.191_joint.vcf"
field_to_visualize = "genotype" # field(s) = "read_depth", *etc.
variant_filter = "pass_only", "list", "significantList_for_proteinstructures.csv"#"range", "chr1:20000000-30000000"
sample_filter = "reorder_columns", "sample_phenotype_matrix.csv", "case_control_status" # "select_columns", "select_column_list.txt",
extension = "pdf"
plot_title = "title_here"




using DataFrames #use CSV.jl ? depwarnings
using CSV
using PlotlyJS
using Rsvg
using Blink
using ViVa

#define functions for backend

function format_reader(vcf,element) #when vcf_matrix is vcf
    format = vcf[1,9]

    s = split(format,":")
    gt_index = find(x -> x == "GT",s)
    dp_index = find(x -> x == "DP",s)
    gq_index = find(x -> x == "GQ",s) #genotype quality - look at annotated pdf to determine how to interpret
    pl_index = find(x -> x == "PL",s)
    mq_index = find(x -> x == "MQ",s)

        if element == "-gt"
            index = gt_index
        elseif element == "-dp"
            index = dp_index
        elseif element == "-gq"
            index = gq_index
        elseif element == "-pl"
            index = pl_index
        elseif element == "-mq"
            index = mq_index
        end

        index = index[1]

    return index

end

function load_vcf(x) #where x is vcf file name to upload 
  readvcf = readlines(x)

  for row = 1:size(readvcf,1)

      if contains(readvcf[row], "#CHROM")
          header_col = row
          global header_col
          header_string = readvcf[row]
      end
  end

  skipstart_number=header_col-1 #this allows vcf to be loaded by readtable, the META lines at beginning of file start with "##" and interfere with readtable function - need readtable vs readdlm for reorder columns
  df_vcf=readtable(x, skipstart=skipstart_number, separator='\t')
    
  vcf=Matrix(df_vcf)

  #load vcf as dataframe twice: once to use for matrix and other to pull header info from
  #df_vcf=readtable(x, skipstart=skipstart_number, separator='\t')

  #2) data cleaning
  ViVa.clean_column1!(vcf)


  for n = 1:size(vcf,1)
      #if typeof(vcf) == "String"
      if vcf[n, 1] != 23 && vcf[n, 1] != 24
      vcf[n, 1] = parse(Int64, vcf[n, 1])
      end
  #end
  end

  #sort rows by chr then chromosome position so are in order of chromosomal architecture
  vcf = sortrows(vcf, by=x->(x[1],x[2]))

  #1) field selection

  #a) FORMAT reader - get index of fields to visualize

    #***once genotype is selected, run format_reader on vcf to get field index
    #global index = format_reader(vcf)
    #return index

return vcf,df_vcf
#return df_vcf #QUESTION for THURSDAY - how to return two variables and define df_vcf when call load_vcf function for use in reorder_columns function
    
end





vcf,df_vcf = load_vcf("variants.filtered.191_joint.vcf")

vcf

println("Upload Successful")
println("Displaying preview of VCF file")
display(vcf)

element = "-gt"
#element = "-gt"

index = format_reader(vcf, element)



#vcf = ViVa.load_sort_phenotype_matrix(phenotype_matrix, key_to_sortby, vcf, df_vcf)
vcf = ViVa.load_sort_phenotype_matrix("sample_phenotype_matrix.csv", "case_control_status", vcf, df_vcf)

vcf = ViVa.select_columns("list_to_include", vcf, df_vcf)

vcf = vcf[(vcf[:,7].== "PASS"),:]

        chr_range = "chromosome_range"

        #create subarray of vcf matching range parameters
        vcf = ViVa.chromosome_range_vcf_filter(chr_range,vcf)


       #load siglist file
        siglist_unsorted=readdlm("variant_list", ',',skipstart=1)

        #replace X with 23 and sort by chr# and position
        ViVa.clean_column1!(siglist_unsorted)
        siglist=sortrows(siglist_unsorted, by=x->(x[1],x[2]))

        #significant variants only filter - subarray of vcf matching lsit of variants of interest
        vcf = ViVa.sig_list_vcf_filter(vcf,siglist)

if element == "-gt"
    
    value_matrix = ViVa.genotype_cell_searcher(vcf, index)

elseif element == "-dp"
    value_matrix = ViVa.dp_cell_searcher(vcf,index)

end

chrlabels=value_matrix[:,1:2]

array_for_plotly = value_matrix[:,10:size(value_matrix,2)]

chr_labeled_array_for_plotly = hcat(chrlabels, array_for_plotly)

array_for_plotly 

writedlm("value_matrix.txt", value_matrix, "\t")

title = "title"

format = "pdf"

if element == "-gt"
    graphic = ViVa.genotype_heatmap2(array_for_plotly, title)
elseif element == "-dp"
    graphic = ViVa.dp_heatmap2(array_for_plotly, title)
end

path = "./"

PlotlyJS.savefig(graphic, "$(path)$title.$(format)")
#JupyterPlot(graphic)

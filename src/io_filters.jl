using GeneticVariation

"""
io_chromosome_range_vcf_filter(x::AbstractString, y::GeneticVariation.VCF.Reader)
create subarray of vcf matching input chromosome range
where x is chromosome range in form chr4:1-3241842
where y is the reader object
"""
function io_chromosome_range_vcf_filter(x::AbstractString, y::GeneticVariation.VCF.Reader)
       a=split(x,":")
       chrwhole=a[1]
       chrnumber=split(chrwhole,"r")
       string_chr=chrnumber[2]
       #chr=parse(string_chr)
       chr=String(string_chr)
       range=a[2]
       splitrange=split(range, "-")
       lower_limit=splitrange[1]
       chr_range_low=parse(lower_limit)
       upper_limit=splitrange[2]
       chr_range_high=parse(upper_limit)
       println(chr_range_high)
       println(chr_range_low)
       println(chr)

       vcf_subarray = Array{Any}(0)

       for record in reader

              if (VCF.chrom(record) == chr) && (chr_range_high > VCF.pos(record) > chr_range_low)
                     #println(record)
                     push!(vcf_subarray,record)
              end
       end

       return vcf_subarray

end

"""
    io_sig_list_vcf_filter(x,y)
Siglist match filter
x is vcf_reader object
y is significant list ordered with substituted chr X/Y for 23/24 from load_siglist()
e.g. function(vcf,siglist)
"""
function io_sig_list_vcf_filter(x::GeneticVariation.VCF.Reader, y)

       vcf_subarray = Array{Any}(0)

       for row= 1:size(y,1)
              dimension = size(y,1)
              #println("dimension is $dimension")

              chr=(y[row,1])
              pos=(y[row,2])

              x = VCF.Reader(open("test_4X_191.vcf", "r"))

              for record in x

                     if typeof(VCF.chrom(record)) == String
                            chr = string(chr)
                            #println(typeof(chr))

#=
                     println(chr)
                     println(VCF.chrom(record))
                     println(typeof(chr))
                     println(typeof(VCF.chrom(record)))
                     println(pos)
                     println(VCF.pos(record))
                     println(typeof(pos))
                     println(typeof(VCF.pos(record)))
=#

                     if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                            #println(record)
                            push!(vcf_subarray,record)
                            #println("match!")
                     end

              else

                     if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                            #println(record)
                            push!(vcf_subarray,record)
                            #println("match!")

              end
              end
              end
       end

       return vcf_subarray
end

"""
io_pass_filter(x::GeneticVariation.VCF.Reader)
create subarray of vcf only including row with FILTER status = PASS
where x is reader object
"""

function io_pass_filter(x::GeneticVariation.VCF.Reader)

       vcf_subarray = Array{Any}(0)

       for record in reader

              #println(VCF.filter(record))

              if VCF.hasfilter(record) && VCF.filter(record) == String["PASS"]
                     #println(record)
                     push!(vcf_subarray,record)

              end
       end
       return vcf_subarray
end

#start of numerical array function
#use @time to see differences between reshape, hcat into matrix vs push! into dataframe
"""
generate_genotype_array(x::Array{Any,1},y::String)
create subarray of vcf only including row with FILTER status = PASS
where x is array of reader objects
where y is genotype field to visualize ex. GT, DP
"""

function generate_genotype_array(x::Array{Any,1},y)

       num_samples = length(VCF.genotype(x[1]))

       df = DataFrame(String, 0, num_samples+2)

       for record in x

              genotype_data_per_variant = VCF.genotype(record, 1:num_samples, y)
              vcf_chrom = VCF.chrom(record)
              vcf_pos = string(VCF.pos(record))
              genotype_data_per_variant = vcat(vcf_pos,genotype_data_per_variant)
              genotype_data_per_variant = vcat(vcf_chrom,genotype_data_per_variant)
              push!(df, genotype_data_per_variant)

       end

       num_array = Matrix(df)
       #ViVa.clean_column1!(num_array)

       return num_array
end

"""
function to create dictionary of values return dict
replace_genotype_with_vals(x::Array)
where x is num_array
"""
function define_geno_dict()
    geno_dict = Dict()

    homo_variant = ["1/1" "1/2" "2/2" "1/3" "2/3" "3/3" "1/4" "2/4" "3/4" "4/4" "1/5" "2/5" "3/5" "4/5" "5/5" "1/6" "2/6" "3/6" "4/6" "5/6" "6/6" "1|1" "1|2" "2|2" "1|3" "2|3" "3|3" "1|4" "2|4" "3|4" "4|4" "1|5" "2|5" "3|5" "4|5" "5|5" "1|6" "2|6" "3|6" "4|6" "5|6" "6|6"]

    hetero_variant = ["0/1" "0/2" "0/3" "0/4" "0/5" "0/6" "1/0" "2/0" "3/0" "4/0" "5/0" "6/0" "0|1" "0|2" "0|3" "0|4" "0|5" "0|6" "1|0" "2|0" "3|0" "4|0" "5|0" "6|0"]

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
function translate_genotype_to_num_array(x::,y::)
function to translate array of genotype to numerical array using dict and genotype array
where x is genotype_array
where y is geno_dict
returns a tuple of num_array for plotting, and chromosome labels for plotting as label bar
"""

function translate_genotype_to_num_array(x,y)

    array_for_plotly = x[:,10:size(x,2)]
    chromosome_labels = x[:,1:2]

    num_array = map(c -> y[c], array_for_plotly)

    return num_array,chromosome_labels
end

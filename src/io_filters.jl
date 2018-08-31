"""
io_chromosome_range_vcf_filter(x::AbstractString, y::GeneticVariation.VCF.Reader)
create subarray of vcf matching input chromosome range
where x is chromosome range in form chr4:1-3241842
where y is the reader object
"""
function io_chromosome_range_vcf_filter(y::GeneticVariation.VCF.Reader, x::AbstractString)
       a=split(x,":")
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

       for record in y

              if (VCF.chrom(record) == chr) && (chr_range_high > VCF.pos(record) > chr_range_low)
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

              chr=(y[row,1])
              pos=(y[row,2])

              x = VCF.Reader(open("test_4X_191.vcf", "r"))

              for record in x

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

"""
io_pass_filter(x::GeneticVariation.VCF.Reader)
create subarray of vcf only including row with FILTER status = PASS
where x is reader object
"""

function io_pass_filter(x::GeneticVariation.VCF.Reader)

vcf_subarray = Array{Any}(0)

for record in x

       if VCF.hasfilter(record) && VCF.filter(record) == String["PASS"]
              push!(vcf_subarray,record)
       end
end

return vcf_subarray
end

"""
generate_genotype_array(x)(x::Array{Any,1},y::String)
create subarray of vcf only including row with FILTER status = PASS
where x is array of reader objects
where y is GT or DP to visualize genotype or read_depth
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

       return num_array
end

"""
       define_geno_dict()
function to create dictionary of values return dict for use in replace_genotype_with_vals())
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
function translate_genotype_to_num_array(x::,y::)
function to translate array of genotype to numerical array using dict and genotype array
where x is genotype_array
where y is geno_dict
returns a tuple of num_array for plotting, and chromosome labels for plotting as label bar
"""
function translate_genotype_to_num_array(x,y)

array_for_plotly = x[:,3:size(x,2)]
chromosome_labels = x[:,1:2]

function translate(c)
       y[c]
end

num_array = map(translate, array_for_plotly)

return num_array,chromosome_labels
end

"""
function translate_readdepth_strings_to_num_array(x::Array{Any,2})
convert read_depth array of strings to int for average calculation
where x is output of generate_genotype_array() for DP option
returns a tuple of num_array type Int for average calculation and plotting, and chromosome labels for plotting as label bar
"""

function translate_readdepth_strings_to_num_array(x::Array{Any,2})

       #clean_column1!(x)

       chromosome_labels = x[:,1:2]

       dp_array_for_plotly = x[:,3:size(x,2)]

       map!(s->replace(s, ".", "0"), dp_array_for_plotly, dp_array_for_plotly)

       dp_array_for_plotly = [parse(Int, i) for i in dp_array_for_plotly]

       for i in dp_array_for_plotly
              if i == "."
                     println(".")
              end
       end

       return dp_array_for_plotly, chromosome_labels
end

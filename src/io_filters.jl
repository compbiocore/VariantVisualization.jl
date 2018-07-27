using GeneticVariation

reader = VCF.Reader(open("test_4X_191.vcf", "r"))

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
              println("dimension is $dimension")

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
                            println("match!")
                     end

              else

                     if (VCF.chrom(record) == chr) && (VCF.pos(record) == pos)
                            #println(record)
                            push!(vcf_subarray,record)
                            println("match!")

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

       if VCF.filter(record) == String["PASS"]
              #println(record)
              push!(vcf_subarray,record)
       end
end

return vcf_subarray
end

n = 0
for record in sub
       n = n+1
       println(VCF.genotype(record, 1:1, "GT"))
       println(n)
end

#start of numerical array function
num_samples = length(VCF.genotype(sub[1]))

numerical_array = Array{Any}(0)

for record in sub

       genotype_data_per_variant = VCF.genotype(record, 1:num_samples, "GT")
       vcat(numerical_array,genotype_data_per_variant)
end

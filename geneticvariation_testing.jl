using GeneticVariation

reader = VCF.Reader(open("test_4X_191.vcf", "r"))

chr_range_low = 1
chr_range_high = 3241842

for record in reader
       if chr_range_high > VCF.pos(record) > chr_range_low
              println(record)
       end
end

#A.) load and clean VCF

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

  #what features to visualize from INFO field - or others - if FORMAT / genotype info isn't included
  #index other data per variant for bars
  #check if format col is included, if so - run function, if not are there cells for each sample? I don't think so - check this.
  global index = format_reader(vcf)

return vcf
return global index

end

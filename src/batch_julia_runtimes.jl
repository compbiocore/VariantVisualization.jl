#!/bin/bash
#=
tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf -l test/test_files/significantList_for_proteinstructures.csv --heatmap genotype`)
println("time to run sig_list heatmap:")
toc()

tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf -r chr1:0-400000000 --heatmap genotype`)
println("time to run chromosome range heatmap:")
toc()

tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf --pass_filter -r chr1:0-400000000 --heatmap genotype`)
println("time to run pass filter and chromosome range heatmap:")
toc()

tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf --pass_filter -x test/test_files/select_column_list.txt --heatmap genotype`)
println("time to run pass filter select 4 patients heatmap:")
toc()
=#
tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf --pass_filter --group_samples test/test_files/sample_phenotype_matrix.csv control,case --heatmap genotype`)
println("time to run pass filter group samples by trait heatmap:")
toc()

tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf -l test/test_files/significantList_for_proteinstructures.csv --heatmap genotype --avg_dp samples`)
println("time to run sig_list heatmap and avg dp scatter plot for samples:")
toc()

tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf -l test/test_files/significantList_for_proteinstructures.csv --heatmap genotype --avg_dp samples`)
println("time to run sig_list heatmap and avg dp scatter plot for variants:")
toc()

tic()
run(`julia viva_cli.jl -f combine_haplo-exo_AC_gatk406.vcf -x test/test_files/select_column_list.txt --heatmap genotype`)
println("time to run select columns with no filters heatmap:")
toc()

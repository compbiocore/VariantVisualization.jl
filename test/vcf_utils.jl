@testset "VCFUtils" begin

#=
#functions to test
clean_column1!,
genotype_cell_searcher_maf_correction *** include this when functional,
genotype_cell_searcher,
dp_cell_searcher,
sig_list_vcf_filter,
chromosome_range_vcf_filter,
load_sort_phenotype_matrix,
reorder_columns,
select_columns
=#

    @testset "clean_column1!" begin
        df = DataFrame(["X" 1 2; "Y" 2 3; 2 4 1])
        clean_column1!(df)
        @test df[1,1] == 23
        @test df[2,1] == 24
        @test df[3,1] == 2
    end

    @testset "genotype_cell_searcher" begin
        df = DataFrame(["X" 1 2; "Y" 2 3; 2 4 1])
        genotype_cell_searcher(x,index) #when x is output subarray of variant selection
        @test df[1,1] == 23
        @test df[2,1] == 24
        @test df[3,1] == 2
    end

end

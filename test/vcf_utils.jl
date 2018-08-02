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

vcf_filename = "test_files/test_4X_191.vcf"

    @testset "clean_column1!" begin

        df = Matrix(["X" 1 2; "Y" 2 3; 2 4 1])
        clean_column1!(df)
        @test df[1,1] == 23
        @test df[2,1] == 24
        @test df[3,1] == 2
    end

    @testset "load_vcf" begin

        vcf_tuple = load_vcf(vcf_filename)
        vcf_df = vcf_tuple[2]
        vcf = vcf_tuple[1]

        @test typeof(vcf_df) == DataFrames.DataFrame
        @test typeof(vcf) == Array{Any,2}
        @test vcf[1,2] == 3076150
        @test vcf_df[1,2] == 3076150
        @test (names(vcf_df))[1] == :_CHROM
        @test (names(vcf_df))[9] == :FORMAT
        @test size(vcf,1) == 1485
        @test size(vcf_df,1) == 1485

        @testset "format_reader" begin

            index_gt = format_reader(vcf,"genotype")
            @test index_gt == 1
            index_dp = format_reader(vcf,"read_depth")
            @test index_dp == 3
        end

        @testset "chromosome_range_vcf_filter(x::AbstractString, vcf::Array)" begin
            chr_range = "chr4:20000000-30000000"
            vcf = chromosome_range_vcf_filter(chr_range,vcf)
            @test vcf[1,1] == 4
            @test vcf[1,2] == 25378769

        end

        vcf = vcf_tuple[1]

        @testset "sig_list_vcf_filter(y::Array,x::Array)" begin

            sig_list_file = "test_files/significantList_for_proteinstructures.csv"
            sig_list = ViVa.load_siglist(sig_list_file)
            vcf = sig_list_vcf_filter(sig_list,vcf)
            @test vcf[1,1] == 4
            @test vcf[1,2] == 3109095

        end

        vcf = vcf_tuple[1]

        #not working - but works in command line - try making sure all items in col_new_order are Symbols
        @testset "load_sort_phenotype_matrix(x::AbstractString, y::AbstractString,vcf::Array,df_vcf::DataFrame)" begin
            pheno_matrix = "test_files/sample_phenotype_matrix.csv"
            vcf = load_sort_phenotype_matrix(pheno_matrix, "case_control_status", vcf, vcf_df)
            @test vcf[1,10] ==  "./.:0,0:0:.:."
            cell_contents = split(vcf[1,10],":")
            @test cell_contents[1] == "./."
            @test vcf[1,2] == 3076150
        end


        @testset "select_columns(x::AbstractString, vcf::Array, df_vcf::DataFrame)" begin

        vcf = select_columns("test_files/select_column_list.txt", vcf, vcf_df)
        @test vcf[1,14] == "./.:0,0:0:.:."
        #println(typeof(vcf[2,10]))
        end

        vcf = vcf_tuple[1]

        @testset "genotype_cell_searcher(x::Array{Any,2}, index::Int64)" begin

        vcf_copy_gt = copy(vcf)
        value_matrix = genotype_cell_searcher(vcf_copy_gt, 1)
        @test value_matrix[1,10] == 0
        @test size(value_matrix) == (1485, 200)
        end

        @testset "dp_cell_searcher(x::Matrix{Any}, index::Int64)" begin

        vcf_copy_dp = copy(vcf)
        value_matrix_dp = dp_cell_searcher(vcf_copy_dp, 3)
        @test value_matrix_dp[1,10] == "0" #does this plot***
        @test size(value_matrix_dp) == (1485, 200)
        end

    end

end


#=***outline for test

#loading and cleaning vcf file
load_vcf(x::AbstractString)

#filter variants
sig_list_vcf_filter(y::Array,x::Array)
chromosome_range_vcf_filter(x::AbstractString, vcf::Array)

#filter samples
load_sort_phenotype_matrix(x::AbstractString, y::AbstractString,vcf::Array,df_vcf::DataFrame)
select_columns(x::AbstractString, vcf::Array, df_vcf::DataFrame)

#field selection and value matrix converter
format_reader(vcf::Array, element::AbstractString)
genotype_cell_searcher(x::Array, index::Int64)
dp_cell_searcher(x::Array{Any,2}, index::Int64)

#analysis
avg_dp_patients(x::Array) #use output array as input for average_dp plotting function


#***plotting tests in separate test set
=#

@testset "VCFUtils" begin
    @testset "clean_column1!" begin
        df = DataFrame(["X" 1 2; "Y" 2 3; 2 4 1])
        clean_column1!(df)
        @test df[1,1] == 23
        @test df[2,1] == 24
        @test df[3,1] == 2
    end



end

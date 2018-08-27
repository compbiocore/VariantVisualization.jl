@testset "PLOTUtils" begin

#=
#functions to test
genotype_heatmap2
dp_heatmap2
=#

@testset "genotype_heatmap2(x,title)" begin

plot = genotype_heatmap2(x,"test_genotype_title")
end

@testset "dp_heatmap2(x,title)" begin

plot = dp_heatmap2(x,"test_read_depth_title")
end

end

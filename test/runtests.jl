module TestViVa

using ViVa.vcf_utils
using ViVa.plot_utils
using DataFrames
using Base.Test


all_tests = [
    ("vcf_utils.jl",   "       Testing: VCF Utils"),
#    ("plot_utils.jl",     "       Testing: Plot Utils")
    ]

println("Running tests:")

for (t, test_string) in all_tests
    println("-----------------------------------------")
    println("-----------------------------------------")
    println(test_string)
    println("-----------------------------------------")
    println("-----------------------------------------")
    include(t)
end

end

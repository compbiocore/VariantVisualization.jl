all_examples = [
    ("literate_src/1_viva.jl",    " Running Example 1")
    ]

println("Running examples:")

for (example, str) in all_examples
    println("-----------------------------------------")
    println("-----------------------------------------")
    println(str)
    println("-----------------------------------------")
    println("-----------------------------------------")

    include(example)
end

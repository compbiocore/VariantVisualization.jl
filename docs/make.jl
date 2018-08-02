using Documenter, ViVa
using Literate

# compile all examples in BioMedQuery/examples/literate_src into markdown and jupyter notebooks for documentation
for (root, dirs, files) in walkdir("examples/literate_src")
    for file in files
        Literate.notebook(joinpath(root,file), joinpath(@__DIR__, "src", "notebooks"))
    end
end

for (root, dirs, files) in walkdir("examples/literate_src")
    for file in files
        Literate.markdown(joinpath(root,file), joinpath(@__DIR__, "src", "examples"))
    end
end

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-material"),
    repo = "github.com/compbiocore/ViVa.jl.git",
    julia  = "0.6",
    osname = "linux"
)

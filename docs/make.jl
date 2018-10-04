using Documenter

makedocs()

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-material"),
    repo = "github.com/compbiocore/ViVa.jl.git",
    julia  = "0.6",
    osname = "linux"
)

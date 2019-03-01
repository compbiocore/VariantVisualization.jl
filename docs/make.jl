using Documenter, VIVA

makedocs(sitename="docs")

deploydocs(
    repo   = "github.com/compbiocore/VIVA.jl.git",
    deps   = Deps.pip("mkdocs", "mkdocs-material", "pygments"),
    make   = () -> run(`mkdocs build`),
    target = "site"
)

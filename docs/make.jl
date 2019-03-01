using Documenter, VIVA

makedocs(sitename="docs")

deploydocs(
    repo   = "github.com/compbiocore/VIVA.jl.git",
    deps   = Deps.pip("mkdocs==1.0.4", "mkdocs-material==4.0.2", "pygments"),
    make   = () -> run(`mkdocs build`),
    target = "site"
)

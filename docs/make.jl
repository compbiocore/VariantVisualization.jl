using Documenter, VIVA

deploydocs(
    repo   = "github.com/compbiocore/VIVA.jl.git",
    deps   = Deps.pip("mkdocs", "pygments", "python-markdown-math"),
    make   = () -> run(`mkdocs build`)
    target = "site"
)

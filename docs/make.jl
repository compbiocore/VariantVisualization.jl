using Documenter, DocumenterMarkdown, VIVA

makedocs(sitename="docs", format = Markdown())

deploydocs(
    repo   = "github.com/compbiocore/VIVA.jl.git",
    deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4", "pygments"),
    make   = () -> run(`mkdocs build`),
    target = "site"
)

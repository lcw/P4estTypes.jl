using P4estTypes
using Documenter

DocMeta.setdocmeta!(P4estTypes, :DocTestSetup, :(using P4estTypes); recursive = true)

makedocs(;
    modules = [P4estTypes],
    authors = "P4estTypes contributors",
    repo = "https://github.com/lcw/P4estTypes.jl/blob/{commit}{path}#{line}",
    sitename = "P4estTypes.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://lcw.github.io/P4estTypes.jl",
        edit_link = "master",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/lcw/P4estTypes.jl",
    devbranch = "master",
    push_preview = true,
)

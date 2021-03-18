using PointSampler
using Documenter

DocMeta.setdocmeta!(PointSampler, :DocTestSetup, :(using PointSampler); recursive=true)

makedocs(;
    modules=[PointSampler],
    authors="Manuel Berkemeier",
    repo="https://github.com/manuelbb-upb/PointSampler.jl/blob/{commit}{path}#{line}",
    sitename="PointSampler.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://manuelbb-upb.github.io/PointSampler.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Thresholded Monte-Carlo Iterator" => "monte_carlo_th.md",
    #    "Docstrings" => "docstrings.md",
    ],
)

deploydocs(;
    repo="github.com/manuelbb-upb/PointSampler.jl",
)
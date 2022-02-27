using Documenter
using GrapheneQFT

push!(LOAD_PATH,"../src/")
makedocs(sitename= "Documentation",
         pages = [
            "Home" => "index.md",
            "Formalism" => "formalism.md",
            "User Guide" => [
               "Installation" => "installation.md",
               "Basics" => "basics.md",
               "Examples" => "examples.md"
            ],
            "API" => "api.md"

         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/rodin-physics/GrapheneQFT.jl.git",
    devbranch = "main"
)

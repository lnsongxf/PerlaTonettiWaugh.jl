## Instructions for Modifying the Package and Compiling Notebooks

## Installation and Use

1. Follow the instructions to [install Julia](https://lectures.quantecon.org/jl/getting_started.html)

2. Install the package, open the Julia REPL (see the documentation above) and then run

```julia
using Pkg
pkg"dev https://github.com/jlperla/PerlaTonettiWaugh.jl.git"
```

The folder will be installed in `~/.julia/dev/PerlaTonettiWaugh` on mac/linux and in somewhere like `C:\Users\USERNAME\.julia\dev\PerlaTonettiWaugh` on Windows

If you are having trouble finding it, then (in the REPL) run
```julia
joinpath(DEPOT_PATH[1], "dev", "PerlaTonettiWaugh")
```

3. Drag that folder into your Git UI (e.g. Github Desktop)

4. Modify the `.jmd` and `.pmd` documents and generate the ipynb with Weave

## Weave Instructions

After installing Weave, run:

```
using IJulia, notebook() # Click yes to the install 
using Conda
Conda.add("nbconvert")
```

To generate the notebooks, run 

```
julia generate.jl
```

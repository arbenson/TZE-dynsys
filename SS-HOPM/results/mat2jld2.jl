# Convert mat to jld2 with v0.6 julia
using MAT
using FileIO, JLD2
function main()
    for file in readdir(".")
        if file[1:length("SS-HOPM")] == "SS-HOPM"
            save(string(split(file, ".")[1], ".jld2"), matread(file))
        end
    end
end

function get_executable_path(package::AbstractString, exec::AbstractString)
    base_path = Base.find_package(package)
    sep = joinpath("src", "$package.jl")
    exec_path = joinpath(split(base_path, sep)[1], exec)

    return exec_path
end

function symlink_user_bin(path_::AbstractString)
    exec = splitdir(path_)[end]
    if Sys.iswindows()
        bin_path = split(path_, "\\$exec")[1]
        @warn bin_path
        run(`setx path "%path%;$bin_path"`)
    else
        ln = "/usr/local/bin/$exec"
        rm(ln, force=true)
        symlink(path_, ln)
        @warn "Created symlink: $ln -> $path_"
    end
end

path_ = get_executable_path("VariantVisualization", "viva")
symlink_user_bin(path_)

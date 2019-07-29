function get_executable_path(package::AbstractString, exec::AbstractString)
    base_path = Base.find_package(package)
    exec_path = joinpath(split(base_path, "src/$package.jl")[1], exec)

    return exec_path
end

function symlink_user_bin(path_::AbstractString)
    exec = split(path_, "/")[end]
    if Sys.iswindows()
        symlink(path_, "C:\\ProgramData\\Bin\\$exec")
    else
        symlink(path_, "/usr/local/bin/$exec")
    end
    @info "Created symlink: /usr/local/bin/viva -> $path_"
end

path_ = get_executable_path("VariantVisualization", "viva")
symlink_user_bin(path_)

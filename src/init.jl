function alias()

    if isfile(joinpath(homedir(), ".zshrc"))
        rcfile = ".zshrc"
    else
        rcfile = ".bashrc"
    end

    dir = split(@__DIR__, "/src")[1]
    exe = joinpath(dir, "viva")

    open(joinpath(homedir(), rcfile), "a+") do io
        write(io, "\n # VIVA Alias \nalias viva=\"$exe\"")
    end
    chmod(exe, 333)
end

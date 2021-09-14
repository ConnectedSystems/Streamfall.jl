using Tar, CodecZlib, ZipFile

const IHACRES_VERSION = v"0.1.0"
const BASE_URL = "https://github.com/ConnectedSystems/ihacres_nim/releases/download/v$(IHACRES_VERSION)/ihacres_nim_"
const HERE = @__DIR__

extract_path = joinpath(HERE, "usr", "lib")

mkpath(extract_path)

if Sys.iswindows()
    target = BASE_URL * "windows.zip"
elseif Sys.islinux()
    target = BASE_URL * "linux.tar.gz"
elseif Sys.isapple()
    target = BASE_URL * "macos.tar.gz"
else
    throw(DomainError("Unsupported platform"))
end

function untar(data::IO, exdir::String)
    if isdir(exdir) == false mkdir(exdir) end
    tar = GzipDecompressorStream(data)
    Tar.extract(tar, exdir; copy_symlinks=true)
    close(tar)
end

fn = download(target, "ihacres_nim")
if endswith(target, "zip")
    dir_name = split(replace(target, ".zip" => ""), "/")[end]
    zarchive = ZipFile.Reader(fn)
    for f in zarchive.files
        if endswith(f.name, "/")
            # skip directories
            continue
        end

        extract_fn = replace(f.name, dir_name * "/" => "")
        file_loc = joinpath(extract_path, extract_fn)
        open(file_loc, "w") do fp
            write(fp, read(f))
        end
    end

    close(zarchive)

elseif endswith(target, "tar.gz")
    dir_name = split(replace(target, ".tar.gz" => ""), "/")[end]
    
    open(fn, "r") do fp
        untar(fp, extract_path)
    end

    mv(joinpath(extract_path, dir_name, "*.*"), joinpath(HERE, "../"), force=true)
end

rm(fn)


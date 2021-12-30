using Tar, CodecZlib, ZipFile, Downloads

const IHACRES_VERSION = v"0.4.1"
const BASE_URL = "https://github.com/ConnectedSystems/ihacres_nim/releases/download/v$(IHACRES_VERSION)/ihacres_nim_"
const HERE = @__DIR__

extract_path = joinpath(HERE, "usr", "lib")

mkpath(extract_path)

if Sys.iswindows()
    target = "windows.zip"
elseif Sys.islinux()
    target = "linux.tar.gz"
elseif Sys.isapple()
    target = "macos.tar.gz"
else
    throw(DomainError("Unsupported platform"))
end

target_url = BASE_URL * target

function untar(data::IO, exdir::String)
    if isdir(exdir) == false mkdir(exdir) end
    tar = GzipDecompressorStream(data)
    Tar.extract(tar, exdir; copy_symlinks=true)
    close(tar)
end

fn = Downloads.download(target_url, "ihacres_nim")
if endswith(target_url, "zip")
    dir_name = split(replace(target_url, ".zip" => ""), "/")[end]
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

elseif endswith(target_url, "tar.gz")
    # dir_name = split(replace(target, ".tar.gz" => ""), "/")[end]

    open(fn, "r") do fp
        untar(fp, abspath(extract_path))
    end

    # mv(joinpath(extract_path, dir_name, "*.*"), joinpath(HERE, "../"), force=true)
end

rm(fn)

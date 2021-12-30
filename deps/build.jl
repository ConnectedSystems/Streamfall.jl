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
    tmp_dir = mktempdir()

    open(fn, "r") do fp
        untar(fp, tmp_dir)
    end

    # Copy files to destination
    # Necessary as Julia's `mv` and `cp` functions
    # require empty directories (`force` option deletes an existing directory)
    for pfn in readdir(tmp_dir)
        open(joinpath(tmp_dir, pfn), "r") do src
            open(joinpath(extract_path, pfn), "w") do dst
                write(dst, read(src))
            end
        end
    end
end

rm(fn)

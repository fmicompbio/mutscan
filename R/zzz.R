.onUnload <- function (libpath) {
    # nocov start
    library.dynam.unload("mutscan", libpath)
    # nocov end
}
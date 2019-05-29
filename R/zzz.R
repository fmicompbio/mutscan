.onUnload <- function (libpath) {
  library.dynam.unload("mutscan", libpath)
}
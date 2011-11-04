.First.lib <- function (lib, pkg) {
    library.dynam("bayesTFR", pkg, lib)
}

.Last.lib <- function (libpath) {
  library.dynam.unload("bayesTFR", libpath)
}

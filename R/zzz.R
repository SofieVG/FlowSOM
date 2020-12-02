.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("Thanks for using FlowSOM. From version 2.1.4 ",
                               "on, the scale \nparameter in the FlowSOM ",
                               "function defaults to FALSE"))
}
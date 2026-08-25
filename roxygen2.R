#
# (1) Make package-development tools available
#

library("devtools")
library("roxygen2")


#
# (2) Create / update documentation files and NAMESPACE
#

# This reads the #' blocks in R/*.R and regenerates man/*.Rd and NAMESPACE.
document()


#
# (3) Build R package and PDF manual
#

build()
build_manual()


#
# (4) Install R package
#

install()


#
# (5) Check R package
#

check()


#
# Notes
#
# For routine documentation updates, document() is the essential command.
# build_manual() requires a working LaTeX installation; it can be skipped if
# only the package tarball and R CMD check are needed.

# Auxiliary functions for morphometrics pipeline
# Authors: João Dornelas & Domingos Cardoso

#_______________________________________________________________________________
# Trait column naming in the analysis (full name vs. code vs. unchanged)

#_______________________________________________________________________________
# Quiet loading of plotting dependencies

# dendextend (also imported by factoextra) and vegan both register a rev()
# method for "hclust" objects, so R prints "Registered S3 method overwritten"
# when the second of them loads. factoextra and dendextend are therefore not
# imported into the moRphoTaxa namespace: they are loaded here, silently, by
# the functions that use them, and called with `pkg::fun()`.
.quiet_load_namespaces <- function(pkgs) {
  for (pkg in pkgs) {
    suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))
  }
  invisible(NULL)
}


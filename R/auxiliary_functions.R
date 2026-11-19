# Auxiliary functions for morphometrics pipeline
# Authors: João Dornelas & Domingos Cardoso


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


#_______________________________________________________________________________
# Function to insert a word after another word in a string ####

insert_after_word <- function(vec, insert_map) {
  for (i in seq_along(insert_map)) {
    after <- names(insert_map)[i]
    insert <- insert_map[[i]]
    if (!insert %in% vec && after %in% vec) {
      vec <- append(vec, insert, after = match(after, vec))
    }
  }
  vec
}

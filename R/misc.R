#' Convert InChIs to SMILES
#'
#' Converts a list of InChIs to SMILES strings.
#'
#' @param inchis strings, list of inchis
#'
#' @return list of smiles strings. The value is "NULL" if the corresponding
#'   InChI could not be converted or was an empty string.
#'
#' @export
inchi2smiles <- function(inchis) {
    smiles <- vector(mode = "character", length = length(inchis))
    for (i in seq(along = inchis)) {
        smi <- system(paste0("obabel -iinchi -:\"", inchis[i], "\" -osmiles"),
                      intern = TRUE, ignore.stderr = TRUE)

        if (length(smi) == 0) {
            # something went wrong when converting the inchi string
            smiles[i] <- "NULL"
        } else {
            smiles[i] <- gsub("[[:space:]]", "", smi)
        }
    }
    return (smiles)
}
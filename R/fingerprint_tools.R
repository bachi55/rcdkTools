#' Calculate molecular fingerprints from inchi
#'
#' This function takes a set of inchis and calculates for those the desired
#' set of fingerprints. For a list of valid fingerprint definitions see
#' \code{\link[rcdk]{get.fingerpint}}.
#'
#' @param inchi list of strings (1 x n_samples)
#' @param fp_type string, fingerprint definition
#' @param fp_mode string, which fingerprint type should be returned:
#' \itemize{
#'   \item "bit": binary fingerprint-vector
#'   \item "raw": (key, value)-pairs
#'   \itemize{
#'     \item key: the ftp definition, e.g. SMARTS string
#'     \item value: the number of occurances of fpt
#'   }
#'   \item "count": counting fingerprint-vector
#' }
#' @param verbose boolean, should the function be verbose
#' @param ... parameters passed to the fingerprint calculation function
#' @return list of \code{\link[fingerprint]{fingerprint-class}} objects for the provided
#'     InChIs. \code{is.null(list[i])} is \code{TRUE}, if the fingerprint for
#'     the corresponding InChI could not be computed.
#'
#' put an example for counting maccs fingerprints here.
#'
#' @export
calculate_fingerprints_from_inchi <- function(
    inchi, fp_type = "maccs", fp_mode = "bit", verbose = FALSE, ...)
{
    smiles <- inchi2smiles(inchi)
    fps <- calculate_fingerprints_from_smiles(smiles, fp_type, fp_mode, verbose, ...)

    stopifnot(names(fps) == smiles)
    names(fps) <- inchi

    return (fps)
}

#' Calculate molecular fingerprints from smiles
#'
#' This function takes a set of smiles and calculates for those the different
#' sets of fingerprints. For a list of valid fingerprint definitions see
#' \code{\link[rcdk]{get.fingerpint}}.
#'
#' @param smiles list of strings (1 x n_samples)
#' @param fp_type string, fingerprint definition
#' @param fp_mode string, which fingerprint type should be returned:
#' \itemize{
#'   \item "bit": binary fingerprint-vector
#'   \item "raw": (key, value)-pairs
#'   \itemize{
#'     \item key: the ftp definition, e.g. SMARTS string
#'     \item value: the number of occurances of fpt
#'   }
#'   \item "count": counting fingerprint-vector
#' }
#' @param verbose boolean, should the function be verbose
#' @param ... parameters passed to the fingerprint calculation function
#' @return list of \code{\link[fingerprint]{fingerprint-class}} objects for the provided
#'     SMILES. \code{is.null(list[i])} is \code{TRUE}, if the fingerprint for
#'     the corresponding SMILES could not be computed.
#'
#' @export
calculate_fingerprints_from_smiles <- function(
    smiles, fp_type = "maccs", fp_mode = "bit", verbose = FALSE, ...)
{
    # Parse all smiles and perform configuration
    tictoc::tic("Parsing and configuration")

    smiles.parsed <- rcdk::parse.smiles (smiles)
    n_mol <- length(smiles)

    # Do not show the progress bar when n_mol <= 1
    verbose <- verbose & n_mol > 1

    if (verbose) {
        pb <- txtProgressBar(0, n_mol - 1, style = 3)
    }

    for (idx in 1:n_mol) {
        # Sometimes the parsing of a smiles failes. For those molecules no
        # descriptors are calculated.
        if (is.null(smiles.parsed[[idx]])) {
            next
        }

        # Do typing for all the molecules
        rcdk::do.typing(smiles.parsed[[idx]])

        if (fp_type == "maccs") {
            # covers also "maccs counting" fingerprints.
            rcdk::do.aromaticity (smiles.parsed[[idx]])
            rcdk::do.isotopes    (smiles.parsed[[idx]]) # ??
        } else if (fp_type == "pubchem") {
            rcdk::do.aromaticity (smiles.parsed[[idx]])
            rcdk::convert.implicit.to.explicit (smiles.parsed[[idx]])
        } else if (fp_type %in% c("circular", "kr", "circular", "signature", "lingo")) {
            # no further configuration needed
        } else if (fp_type == "substructure") {
            # Unfortunately, we need do not know here, whether the user has
            # provided here some substructures that would require some kind of
            # modifications.
        } else if (fp_type == "estate") {
            rcdk::do.aromaticity (smiles.parsed[[idx]])
            # rJava::.jcall(smiles.parsed[[idx]], "V", "addImplicitHydrogens", smiles.parsed[[idx]])
        } else if (fp_type %in% c("standard", "extended", "graph", "hybridization")) {
            # All based on the 'Fingerprinter' class
            rcdk::convert.implicit.to.explicit (smiles.parsed[[idx]])
            # aromaticity detection is done in the fingerprint class (CDK)
        } else if (fp_type == "shortestpath") {
            rcdk::convert.implicit.to.explicit (smiles.parsed[[idx]])
            # aromaticity detection is done in the fingerprint class (CDK)
        } else {
            stop("Unsupported fingerprint type: ", fp_type)
        }

        if (verbose) {
            setTxtProgressBar(pb, idx)
        }
    }

    if (verbose) {
        close(pb)
    }

    tictoc::toc(log = TRUE, quiet = ! verbose)

    # Calculate all the desired fingerprints
    if (verbose) {
        tictoc::tic(sprintf("\nCalculate the '%s' fingerprints", fp_type))
        pb <- txtProgressBar(0, n_mol - 1, style = 3)
    }

    fps <- lapply(smiles.parsed,
        function(x) {
            if (is.null(x)) {
                # SMILES could not be parsed
                NULL
            } else {
                if (verbose){
                    setTxtProgressBar(pb, which(sapply (smiles.parsed, '==', x)))
                }

                if ((fp_type == "maccs") & (fp_mode == "count")) {
                    # This is a bit hacky, as it ignores the parameters passed by
                    # using '...'. However, at this point, the 'substructure' fps
                    # do anyway not accept any other parameters.
                    rcdk::get.fingerprint(x, type = "substructure", fp.mode = fp_mode,
                                          substructure.pattern = count_maccs_pattern)
                } else {
                    rcdk::get.fingerprint(x, type = fp_type, fp.mode = fp_mode, ...)
                }

            }
        })

    if (verbose) {
        tictoc::toc(log = TRUE, quiet = ! verbose)
        close(pb)
    }

    return (fps)
}

#' Convert fingerprint list into a JSON string and write it to a file.
#'
#' List of \code{\link[fingerprint]{fingerprint}} objects is converted to a JSON
#' string. Subsequently this string is stored in a file using
#' \code{\link[jsonlite]{write_json}}.
#'
#' @param fps list of \code{\link[fingerprint]{fingerprint}} objects
#'   (1 x n_samples)
#' @param path string, path to the output file
#' @param exclude_zero_fp logical, indicating whether fingerprints that have
#'   value equal 0 should be excluded from the JSON string.
#' @param ... parameters passed to down to \code{\link[jsonlite]{toJSON}}.
#'
#' @export
write_fingerprint_to_json_file <- function (fps, path, exclude_zero_fp = FALSE,
                                            ...)
{
    fps_out <- lapply(fps, function(fp_i) {
        if (class(fp_i) == "featvec")
        {
            features <- sapply(fp_i@features, fingerprint::count)
            names(features) <- sapply(fp_i@features, fingerprint::feature)
        }
        else if (class(fp_i) == "fingerprint")
        {
            features <- rep(0, fp_i@nbit)
            features[fp_i@bits] <- 1
            names(features) <- 1:fp_i@nbit
        }
        else stop("Unsupported fingerprint class: ", class(fp_i))

        if (exclude_zero_fp) {
            features <- features[features > 0]
        }

        lapply(features, function (x) x)
    })

    jsonlite::write_json(fps_out, path, ...)
}


#' Convert fingerprint list into a matrix and write it to a csv-file.
#'
#' List of \code{\link[fingerprint]{fingerprint}} objects is converted to a
#' matrix using \code{\link{fingerprint_to_matrix}} and subsequently
#' written into a csv-file.
#'
#' @param fps list of \code{\link[fingerprint]{fingerprint}} objects
#'   (1 x n_samples)
#' @param path string, path to the output file
#'
#' @export
write_fingerprint_to_csv_file <- function (fps, path) {
    if (class(fps) != "matrix") {
        fps_matrix <- fingerprints_to_matrix(fps)
    } else {
        fps_matrix <- fps
    }
    write.table(fps_matrix, path, col.names = FALSE, row.names = TRUE,
                quote = TRUE, sep = ",")
}

#' Construct a fingerprint matrix
#'
#' List of \code{\link[fingerprint]{fingerprint}} objects to matrix.
#'
#' @param fps list of \code{\link[fingerprint]{fingerprint}} objects
#'   (1 x n_samples)
#' @return fingerprint matrix (n_samples x n_fingerprints)
#'
#' @export
fingerprints_to_matrix <- function (fps)
{
    fps_matrix <- fp.to.matrix2(fps)

    # get InChIs or SMILES associated with each fingerprint.
    rownames(fps_matrix) <- names(fps)

    return (fps_matrix)
}

fp.to.matrix2 <- function (fps) {
    # get_value_from_feature <- function (feature) { attributes(feature)$count }
    # get_key_from_feature <- function (feature) { attributes(feature)$feature }

    # Determine fingerprint class
    is_binary_fp <- NULL
    for (fp in fps) {
        if (! is.null(fp)) {
            is_binary_fp <- class(fp) == "fingerprint"
            break
        }
    }
    if (is.null(is_binary_fp)) {
        stop("All provided fingerprints are NULL.")
    }

    # Determine the fingerprint dimension
    fp_dims <- get_fingerprint_dimensions(fps, check_all = TRUE, silent = TRUE)
    if (length(fp_dims) > 1) {
        stop("Fingerprint vectors do have different dimension. Cannot create a
             single matrix")
    }

    m <- matrix (0, nrow = length(fps), ncol = fp_dims)
    for (fp_i in seq(1, length(fps))) {
        fp <- fps[[fp_i]]

        # for NULL fps we insert NA for each fp-dimension fp_i
        if (is.null(fp)) {
            m[fp_i, ] <- NA
            next
        }

        if (is_binary_fp) {
            m[fp_i, fp@bits] <- 1
        } else {
            m[fp_i, ] <- sapply(fp@features, fingerprint::count)
        }
    }

    return (m)
}

#' Determine the fingerprint dimension
#'
#' Given a list of fingerprints, e.g. the output of
#' \link{calculate_fingerprints_from_inchi}, the dimension of the fingerprints
#' is determined.
#'
#' @param fps list of \code{\link[fingerprint]{fingerprint}} object
#' @param check_all logical, if \code{FALSE} only the first not null fingerprint
#'   is used to determine the fingerprint dimension.
#' @param silent logical, if \code{TRUE} warnings are supressed.
#'
#' @return integer, dimension of the fingerprints or list of integers, if several
#'   different dimensions where found.
#'
#' @export
get_fingerprint_dimensions <- function (fps, check_all = FALSE, silent = FALSE) {
    fp_dims <- numeric()
    for (fp in fps) {
        if (! is.null(fp)) {
            if (class(fp) == "fingerprint") {
                fp_dim <- attributes(fp)$nbit
            } else {
                fp_dim <- length(attributes(fp)$features)
            }

            if (is.null(fp_dim)) {
                stop("Fingerprint does have no 'nbit' or 'features' attribute.")
            }

            if (check_all) {
                fp_dims <- c(fp_dims, fp_dim)
            } else {
                fp_dims <- fp_dim
                break
            }
        }
    }

    # If all fingerprints are NULL, we cannot determine the dimension.
    if (! silent && length(fp_dims) == 0) {
        warning("All fingerprints are 'NULL'. No dimension can be determined.")
        fp_dims <- 0
    }

    fp_dims <- unique(fp_dims)

    if (! silent && length(fp_dims) > 1) {
        warning("Fingerprints do have different dimensions.")
    }

    return(fp_dims)
}

#' Calculate mask to exclude molecular fingerpints given a set of criterias
#'
#' Given a binary fingerprint matrix a mask is calculated to exclude
#' fingerprints:
#' \itemize{
#'   \item with the same value across the dataset
#'   \item with low variance, e.g. fp_i = 1 in 90\% of the data
#'   \item which are redunant, e.g. fp_i = 1 <=> fp_j = 1
#' }
#' where fp_i is the i'th fingerprint definition.
#'
#' @param fps binary matrix, shape (n_samples x n_fingerprints), e.g. the output
#'   of \code{fingerprints_to_matrix}.
#' @param remove_single_value binary, exclude fps with always the same value
#' @param remove_redundant binary, exclude redundant fps
#' @param remove_low_variance binary, exlcude fps with low variance
#' @param low_variance_tshd scalar, threshold for low variance (default = 0.9,
#'   i.e. fps where >=90\% of the examples have the same value are removed)
#'
#' @return binary vector (1 x n_fingerpints):
#' \itemize{
#'   \item \code{TRUE}: keep fingerprint
#'   \item \code{FALSE}: exclude fingerprint
#' }
#'
#' @export
get_fingerprint_mask <- function (
    fps,
    remove_single_value = TRUE,
    remove_redundant    = TRUE,
    remove_low_variance = TRUE,
    low_variance_tshd   = 0.90)
{
    stopifnot(all(unique(as.vector(fps)) %in% c(TRUE, FALSE, 1, 0)))

    n_samples <- nrow(fps)
    n_fps <- ncol(fps)

    # Find the columns in which all the fingerprints are either 0 (FALSE) or 1
    # (TRUE).
    if (remove_single_value) {
        is_single_value <- apply(fps, MARGIN = 2,
                                 FUN = function(x) all(x == 1) || all(x == 0))
    } else {
        is_single_value <- rep(FALSE, n_fps)
    }

    # Find columns with low variance of the fingerprints.
    if (remove_low_variance) {
        is_low_variance <- apply(fps, MARGIN = 2, FUN = function(x) {
            any(c(sum(x == 1), sum(x == 0)) / n_fps >= low_variance_tshd)
        })
    } else {
        is_low_variance <- rep(FALSE, n_fps)
    }

    # Find all the redundant fingerprints:
    #   fp_i = 0 <=> fp_j = 0
    #   fp_i = 1 <=> fp_j = 1
    if (remove_redundant) {
        fps[fps == 0] <- -1

        fps_cor <- (t(fps) %*% fps) / n_samples

        # Build up an undirected graph with edges between all redundant
        # components
        fps_graph <- igraph::graph_from_adjacency_matrix(fps_cor == 1,
                                                         mode = "undirected")
        # Find the connected components in this graph. For each component only
        # one fingerprint definition needs to be kept. We keep the one with the
        # lowest fingerprint index.
        fps_comp <- igraph::components(fps_graph)
        n_comp <- fps_comp$no
        cols_not_redundant <- sapply(1:n_comp, FUN = function (idx, fps_comp) {
                which (fps_comp$membership == idx)[1]
            }, fps_comp)

        is_redundant <- rep (TRUE, n_fps)
        is_redundant[cols_not_redundant] <- FALSE
    } else {
        is_redundant <- rep (FALSE, n_fps)
    }

    return (! (is_single_value | is_low_variance | is_redundant))
}


#' Calculate mask to exclude molecular counting fingerpints
#'
#' @param fps scalar matrix, shape (n_samples x n_fingerprints), e.g. the output
#'   of \code{fingerprints_to_matrix}.
#' @param remove_single_value binary, exclude fps with always the same value
#'
#' @return binary vector (1 x n_fingerpints):
#' \itemize{
#'   \item \code{TRUE}: keep fingerprint
#'   \item \code{FALSE}: exclude fingerprint
#' }
#'
#' @export
get_count_fingerprint_mask <- function (
    count_fps,
    remove_single_value = TRUE)
{
    n_fps <- ncol(count_fps)

    if (remove_single_value) {
        is_single_value <- apply (count_fps, MARGIN = 2,
                                  FUN = function (x) length(unique(x)) == 1)
    } else {
        is_single_value <- rep(FALSE, n_fps)
    }

    return (! is_single_value)
}
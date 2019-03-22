#' Calculate molecular descriptors from InChI
#'
#' This function takes a set of InChI and calculates for those the different
#' desired molecular descriptors provided by the CDK library through rcdk.
#'
#' @param inchi list of strings (1 x n_samples)
#' @param desc_name string, name of the descriptor to calculate
#' @param skip_3D_desc boolean, indicating whether descriptors that require
#'     3D coordinates should be skipped. Return value of the function will be
#'     \code{NULL} in this case.
#' @param verbose boolean, should the function be verbose
#' @return named list
#' \itemize{
#'   \item \code{desc_vals}, named list with the descriptor vectors for each InChI
#'   \item \code{d_type}, string data type of the descriptor
#'   \item \code{desc_params}, named list of descriptor parameters used
#' }
#' @export
#'
#' @export
calculate_descriptor_from_inchi <- function(
    smiles, desc_name, skip_3D_desc=FALSE, verbose=FALSE)
{
    smiles <- inchi2smiles(inchi)
    descs <- calculate_descriptor_from_smiles(smiles, desc_name, skip_3D_desc, verbose)

    stopifnot(names(descs) == smiles)
    names(descs) <- inchi

    return (descs)
}

#' Calculate molecular descriptors from SMILES
#'
#' This function takes a set of SMILES and calculates for those the different
#' desired molecular descriptors provided by the CDK library through rcdk.
#'
#' @param smiles list of strings (1 x n_samples)
#' @param desc_name string, name of the descriptor to calculate
#' @param skip_3D_desc boolean, indicating whether descriptors that require
#'     3D coordinates should be skipped. Return value of the function will be
#'     \code{NULL} in this case.
#' @param verbose boolean, should the function be verbose
#' @return named list
#' \itemize{
#'   \item \code{desc_vals}, named list with the descriptor vectors for each SMILES
#'   \item \code{d_type}, string data type of the descriptor
#'   \item \code{desc_params}, named list of descriptor parameters used
#' }
#' @export
calculate_descriptor_from_smiles <- function(
    smiles, desc_name, skip_3D_desc=TRUE, verbose=FALSE)
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

    # Get meta data for the requested descriptor
    desc_meta <- get_descriptor_meta_data(desc_name)
    if (skip_3D_desc & desc_meta$is_3D) {
        warning("3D descriptors will be skipped.")
        return(NULL)
    }

    for (idx in 1:n_mol) {
        # Sometimes the parsing of a smiles failes. For those molecules no
        # descriptors are calculated.
        if (is.null(smiles.parsed[[idx]])) {
            next
        }

        # Apply the required configuration
        for (fun in desc_meta$conf_func){
            fun(smiles.parsed[[idx]])
        }

        if (verbose) {
            setTxtProgressBar(pb, idx)
        }
    }

    if (verbose) {
        close(pb)
    }

    tictoc::toc(log = TRUE, quiet = ! verbose)

    # Calculate all the desired descriptors
    if (verbose) {
        tictoc::tic(sprintf("\nCalculate the '%s' descriptor", desc_name))
        pb <- txtProgressBar(0, n_mol - 1, style = 3)
    }

    desc_vals <- lapply(smiles.parsed,
    function(x) {
        if (is.null(x)) {
            # SMILES could not be parsed
            NULL
        } else {
            if (verbose){
                setTxtProgressBar(pb, which(sapply (smiles.parsed, '==', x)))
            }

            # if (desc_meta$clone_atom_container) {
            #     x_copy = rJava::.jcall(
            #         x, "Lorg/openscience/cdk/interfaces/IAtomContainer;", "clone")
            # } else {
            #     x_copy = x
            # }

            loc_desc_vals <- rcdk::eval.desc(x, desc_name,
                                             desc.params = desc_meta$desc_params)

            stopifnot(length(loc_desc_vals) == desc_meta$d_length)

            return(as.numeric(loc_desc_vals))
        }
    })

    if (verbose) {
        tictoc::toc(log = TRUE, quiet = ! verbose)
        close(pb)
    }

    return (list(desc_vals=desc_vals, d_type=desc_meta$d_type,
                 desc_params=desc_meta$desc_params, d_length=desc_meta$d_length))
}

#' Metadata for CDK descriptors
#'
#' This function returns the metadata, e.g. needs 3D, is sclar, etc., for the
#' specified descriptor. This information is mainly extracted from the CDK java
#' code.
#'
#' @param desc_name string, rcdk-like descriptor name
#'
#' @return named list containing the descriptor metadata:
#' \itemize{
#'   \item \code{is_3D}, boolean indicating whether the descriptor requires 3D coordinates
#'   \item \code{conf_func}, list of rcdk molecule-configuration functions to be
#'     applied before descriptor calculation
#'   \item \code{clone_atom_container}, boolean indicating, whether the atom
#'     container should be cloned before the descriptor is calculated.
#'   \item \code{desc_params}, vector of names lists that can be passed to the descriptor
#'   \item \code{d_type}, string either float or int(eger)
#'   \item \code{d_length}, length of the descriptor vector
#' }
get_descriptor_meta_data <- function(desc_name) {
    meta_data = list(is_3D=FALSE,
                     conf_func=c(),
                     clone_atom_container=FALSE,
                     desc_params=list(),
                     d_type="",
                     d_length=-1)

    if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor"){
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor"){
        # list of integers (counts)
        meta_data$d_type <- "int"
        meta_data$d_length <- 11
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor"){
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor"){
        # scalar integer (sum of squared counts)
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor"){
        # scalar float
        meta_data$conf_func <- c(rcdk::do.typing,
                                 rcdk::do.aromaticity,
                                 rcdk::convert.implicit.to.explicit)
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor"){
        # list of integer (sum of distances)
        meta_data$d_type <- "int"
        meta_data$d_length <- 2
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor"){
        # list of floats
        meta_data$is_3D <- TRUE
        meta_data$d_type <- "float"
        meta_data$d_length <- 17
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 5
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor") {
        # list of floats
        meta_data$is_3D <- TRUE
        meta_data$d_type <- "float"
        meta_data$d_length <- 2
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor") {
        meta_data$is_3D <- TRUE
        meta_data$d_type <- "float"
        meta_data$d_length <- 7
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 19
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$clone_atom_container <- TRUE
        meta_data$desc_params <- list(checkRingSystem=TRUE)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor") {
        # list of floats
        meta_data$is_3D <- TRUE
        meta_data$d_type <- "float"
        meta_data$d_length <- 2
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$clone_atom_container <- TRUE
        meta_data$desc_params <- list(checkRingSystem=TRUE, checkAromaticity=TRUE)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor") {
        # list of counts
        meta_data$d_type <- "int"
        meta_data$d_length <- 79
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 3
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor") {
        # list of floats
        meta_data$is_3D <- TRUE
        meta_data$d_type <- "float"
        meta_data$d_length <- 9
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor") {
        # list of floats
        meta_data$is_3D <- TRUE
        meta_data$d_type <- "float"
        meta_data$d_length <- 29
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 16
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 6
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 8
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 10
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor") {
        # list of "types" (integer)
        meta_data$d_type <- "int"
        meta_data$d_length <- 9
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1

        # TODO: How to handle the different parameter sets?
        # meta_data$d_length <- 4
        # meta_data$desc_params <- c(list(order="s"), list(order="d"),
        #                            list(order="t"), list(order="q"))
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 6
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 5
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 5
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 5
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor") {
        # scalar float
        meta_data$d_type <- "float"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::convert.implicit.to.explicit)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor") {
        # list of floats
        meta_data$d_type <- "float"
        meta_data$d_length <- 3
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity,
                                 rcdk::convert.implicit.to.explicit)
    } else if (desc_name == "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor") {
        # scalar count
        meta_data$d_type <- "int"
        meta_data$d_length <- 1
        meta_data$conf_func <- c(rcdk::do.typing, rcdk::do.aromaticity)
    } else {
        stop(paste("Invalid descriptor name:", desc_name))
    }

    return(meta_data)
}

#' Construct a descriptor matrix from \code{calculate_descriptor_from_*} function
#'
#' Given a list of descriptor values, i.e., the values of a _single_ descriptor
#' for a _set_ of molecules, this function will return a matrix object containing
#' all values.
#'
#' @param desc Output of \code{calculate_descriptor_from_*} functions, i.e. the
#'    descriptor values for a set of (n_samples, ) molecules.
#'
#' @return fingerprint matrix (n_samples x n_desc_vals)
#'
#' @export
descriptors_to_matrix <- function (desc)
{
    # Number of molecules: n_samples
    n_samples <- length(desc$desc_vals)

    # Number of values associated with the descriptor: n_desc_vals
    n_desc_vals <- desc$d_length

    # Create output matrix
    desc_matrix <- matrix(nrow=n_samples, ncol=n_desc_vals)
    rownames(desc_matrix) <- names(desc$desc_vals)

    for(i_sample in seq(along=desc$desc_vals)) {
        if(is.null(desc$desc_vals[[i_sample]])) {
            next
        }
        desc_matrix[i_sample, ] <- desc$desc_vals[[i_sample]]
    }

    return (desc_matrix)
}
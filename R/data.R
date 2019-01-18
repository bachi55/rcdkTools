#' MACCS SMARTS pattern that can be counted
#'
#' List of SMARTS pattern extracted from the MACCS fingerprint definition that
#' are countable. Those patterns are used of the "maccs" fingerprints with the
#' mode "count" should be calcualted, e.g. using
#' \code{\link{calculate_fingerprints_from_inchi}}.
#'
#' @format list of strings, SMARTS pattern
"count_maccs_pattern"

#' MACCS SMARTS pattern that represent binary properties.
#'
#' List of SMARTS pattern extracted from the MACCS fingerprint definition that
#' are binary by nature. Those patterns can be calculated using 'substructure'
#' fingerprints with specified \code{substructure.pattern}. Please check
#' \code{\link[rcdk]{get.fingerprint}} for further details.
#'
#' @format list of strings, SMARTS pattern
"binary_maccs_pattern"

#' OpenBabel FP3 SMARTS pattern
#'
#' List of SMARTS pattern provided by OpenBabel [1]. The pattern capture
#' functional groups of molecules using 55 SMARTS:
#'
#' 1 cation
#' 2 anion
#' 3 aldehyde or ketone
#' 4 aldehyde
#' 5 ketone
#' 6 thioaldehyde or thioketone
#' 7 thioaldehyde
#' 8 thioketone
#' 9 imine
#' 10 hydrazone
#' 11 semicarbazone
#' 12 thiosemicarbazone
#' 13 oxime
#' 14 oxime ether
#' 15 ketene
#' 16 keten acetyl derivative***
#' 17 carbonyl hydrate
#' 18 hemiacetal
#' 19 acetal
#' 20 hemiaminal
#' 21 aminal
#' 22 thiohemiaminal
#' 23 thioacetal
#' 24 enamine
#' 25 enol
#' 26 enol ether
#' 27 hydroxy compound
#' 28 alcohol
#' 29 primary alcohol
#' 30 secondary alcohol
#' 31 tertiary alcohol
#' 32 1,2-diol
#' 33 1,2-aminoalcohol
#' 34 phenol
#' 35 1,2-diphenol
#' 36 enediol
#' 37 ether
#' 38 dialkyl ether
#' 39 alkylaryl ether
#' 40 diaryl ether
#' 41 thioether
#' 42 disulfide
#' 43 peroxide
#' 44 hydroperoxide
#' 45 aryl
#' 46 heteroatom
#' 47 HBA
#' 48 HBD
#' 49 Atom in a ring
#' 50 carboxylic acid
#' 51 ester
#' 52 nitro
#' 53 nitrile
#' 54 aniline
#' 55 urea
#'
#' The patters originate from a list of roughly 200 functional groups [2]. To
#' visualize the SMARTS pattern you can use the following webservice [3].
#'
#' @references
#' \itemize{
#'   \item [1] \url{http://openbabel.org/wiki/Tutorial:Fingerprints}
#'   \item [2] \url{http://merian.pch.univie.ac.at/~nhaider/cheminf/fgtable.pdf}
#'   \item [3] \url{https://smarts.plus/}
#' }
#'
#' @format list of strings, SMARTS pattern
"fp3_obabel_pattern"

#' OpenBabel logp SMARTS pattern
#'
#' List of atom type SMARTS classification patterns extracted from [1].
#'
#' @references
#' \itemize{
#'   \item [1] \url{https://doi.org/10.1021/ci990307l}
#' }
#'
#' @format list of strings, SMARTS pattern
"logp_obabel_pattern"
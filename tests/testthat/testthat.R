library(rcdkTools)

context("Convert molecular representations")
test_that("inchis_are_correctly_converted", {
    inchis <- c(
        "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3/t4-/m1/s1",
        "InChI=1S/C4H9Cl/c1-3-4ydfdsfsd",
        "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3",
        "",
        "InChI=1S/C20H19F3N2O4/c1-13(14-8-6-9-16(11-14)20(21,22)23)24-29-12-15-7-4-5-10-17(15)18(25-28-3)19(26)27-2/h4-11H,12H2,1-3H3/b24-13+,25-18-",
        "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3",
        "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)")

    smiles <- inchi2smiles(inchis)

    expect_equal(length(smiles), length(inchis))
    expect_equal(smiles[1], "CC[C@@H](C)Cl")
    expect_equal(smiles[2], "NULL")
    expect_equal(smiles[3], "CCC(C)Cl")
    expect_equal(smiles[4], "NULL")
    expect_equal(smiles[5], "C/C(=N\\OCc1ccccc1/C(=N/OC)/C(=O)OC)/c1cccc(c1)C(F)(F)F")
    expect_equal(smiles[6], "CCC(C)Cl")
    expect_equal(smiles[7], "CC(C(=O)O)N")
})

context("Binary fingerprint calculation from SMILES")
test_that("fingerprints are calculated", {
    smiles <- c("CCC(C)Cl", "CC(C(=O)O)N", "C/C(=N\\OCc1ccccc1/C(=N/OC)/C(=O)OC)/c1cccc(c1)C(F)(F)F")
    fps <- calculate_fingerprints_from_smiles(smiles)

    expect_equal(length(fps), length(smiles))
    expect_equal(names(fps), smiles)
    for (fp in fps) {
        expect_false(is.null(fp))
    }

    fps_matrix <- fingerprint::fp.to.matrix(fps)
    expect_equal(ncol(fps_matrix), 166)
})

test_that("parameters are passed to fingerprint calculation", {
    smiles <- c("CCC(C)Cl", "NC(C)C(=O)O", "F/C=C/F")
    fps <- calculate_fingerprints_from_smiles(smiles, fp_type = "circular")

    expect_equal(length(fps), length(smiles))
    expect_equal(names(fps), smiles)
    for (fp in fps) {
        expect_false(is.null(fp))
    }

    fps_matrix <- fingerprint::fp.to.matrix(fps)
    expect_equal(ncol(fps_matrix), 1024)

    for (size in c(128, 256, 512, 1024)) {
        fps <- calculate_fingerprints_from_smiles(smiles, fp_type = "standard",
                                                  size = size)

        expect_equal(length(fps), length(smiles))
        expect_equal(names(fps), smiles)
        for (fp in fps) {
            expect_false(is.null(fp))
        }

        fps_matrix <- fingerprint::fp.to.matrix(fps)
        expect_equal(ncol(fps_matrix), size)
    }
})

test_that("non parseable smiles are handled", {
    smiles <- c("CCC(C)Cl", "NC(C)C(=Oasd)O", "F/C=C/Fasdasd")
    fps <- calculate_fingerprints_from_smiles(smiles)

    expect_equal(length(fps), length(smiles))
    expect_equal(names(fps), smiles)

    expect_false(is.null(fps[[1]]))
    expect_true(is.null(fps[[2]]))
    expect_true(is.null(fps[[3]]))
})

context("Binary fingerprint calculation from InChIs")
test_that("fingerprints are calculated", {
    inchis <- c(
        "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3/t4-/m1/s1",
        "InChI=1S/C4H9Cl/c1-3-4ydfdsfsd",
        "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3",
        "",
        "InChI=1S/C20H19F3N2O4/c1-13(14-8-6-9-16(11-14)20(21,22)23)24-29-12-15-7-4-5-10-17(15)18(25-28-3)19(26)27-2/h4-11H,12H2,1-3H3/b24-13+,25-18-",
        "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3",
        "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)")

    smiles <- inchi2smiles(inchis)

    fps <- calculate_fingerprints_from_inchi(inchis)
    fps_from_smiles <- calculate_fingerprints_from_smiles(smiles)


    expect_equal(names(fps), inchis)

    # Except for the names, the calculated fingerprints must be the same.
    names(fps) <- NULL
    names(fps_from_smiles) <- NULL
    expect_equal(fps, fps_from_smiles)

    expect_false(is.null(fps[[1]]))
    expect_true(is.null(fps[[2]]))
    expect_false(is.null(fps[[3]]))
    expect_true(is.null(fps[[4]]))
    expect_false(is.null(fps[[5]]))
    expect_false(is.null(fps[[6]]))
    expect_false(is.null(fps[[7]]))
})

test_that("parameters are passed to fingerprint calculation", {
    inchis <- c("InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3/t4-/m1/s1",
                "InChI=1S/C4H9Cl/c1-3-4(2)5/h4H,3H2,1-2H3")
    fps <- calculate_fingerprints_from_inchi(inchis, fp_type = "circular")

    expect_equal(length(fps), length(inchis))
    expect_equal(names(fps), inchis)
    for (fp in fps) {
        expect_false(is.null(fp))
    }

    fps_matrix <- fingerprint::fp.to.matrix(fps)
    expect_equal(ncol(fps_matrix), 1024)

    for (size in c(128, 256, 512, 1024)) {
        fps <- calculate_fingerprints_from_inchi(inchis, fp_type = "standard",
                                                  size = size)

        expect_equal(length(fps), length(inchis))
        expect_equal(names(fps), inchis)
        for (fp in fps) {
            expect_false(is.null(fp))
        }

        fps_matrix <- fingerprint::fp.to.matrix(fps)
        expect_equal(ncol(fps_matrix), size)
    }
})

context("Create binary fingerprint matrix")
test_that("matrix is equal to the result of fingerprint::fp.to.matrix", {
    smiles <- c("CCC(C)Cl", "NC(C)C(=O)O", "NULL", "F/C=C/F")
    fps <- calculate_fingerprints_from_smiles(smiles, fp_type = "circular")

    fps_matrix <- fingerprints_to_matrix(fps)

    expect_equal(rownames(fps_matrix), smiles)
    expect_equal(nrow(fps_matrix), 4)
    expect_equal(ncol(fps_matrix), 1024)

    expect_false(any(is.na(fps_matrix[-3, ])))
    expect_true(all(is.na(fps_matrix[3, ])))

    expect_equal(which(fps_matrix[1,]==1), attributes(fps[[1]])$bits)
    expect_equal(which(fps_matrix[2,]==1), attributes(fps[[2]])$bits)
    expect_equal(which(fps_matrix[4,]==1), attributes(fps[[4]])$bits)
})

test_that("some obscure cases are caught correctly", {
    fps <- list(new("fingerprint", nbit=6, bits=c(1,2,5,6)),
                new("fingerprint", nbit=8, bits=c(1,4,5,6)),
                new("fingerprint", nbit=6, bits=c(2,3,4,5,6)))
    expect_error(fingerprints_to_matrix(fps))

    fps <- list(NULL, NULL)
    expect_error(fingerprints_to_matrix(fps))
})

context("Counting fingerprint mask")
test_that("is_correct", {
    fps <- matrix (c(
        1, 6, 3, 0,
        2, 0, 3, 1,
        0, 0, 3, 1
    ), nrow = 3, byrow = TRUE)

    mask <- get_count_fingerprint_mask (fps, remove_single_value = FALSE,
                                        remove_low_abundant = FALSE)
    expect_equal(mask, c(T, T, T, T))

    mask <- get_count_fingerprint_mask (fps, remove_single_value = TRUE,
                                        remove_low_abundant = FALSE)
    expect_equal(mask, c(T, T, F, T))
})

test_that("low_abundance_detection_is_correct", {
    fps <- matrix (c(
        1, 6, 3, 0,
        2, 0, 3, 0,
        0, 0, 3, 0,
        0, 0, 3, 1
    ), nrow = 4, byrow = TRUE)

    mask <- get_count_fingerprint_mask (fps, remove_low_abundant = FALSE,
                                        remove_single_value = FALSE)
    expect_equal(mask, c(T, T, T, T))

    mask <- get_count_fingerprint_mask (fps, remove_low_abundant = TRUE,
                                        remove_single_value = FALSE,
                                        low_abundance_thsd = 0.30)
    expect_equal(mask, c(T, F, T, F))

    mask <- get_count_fingerprint_mask (fps, remove_low_abundant = TRUE,
                                        remove_single_value = FALSE,
                                        low_abundance_thsd = 0.51)
    expect_equal(mask, c(F, F, T, F))

    mask <- get_count_fingerprint_mask (fps, remove_low_abundant = TRUE,
                                        remove_single_value = FALSE)
    expect_equal(mask, c(T, T, T, T))

    mask <- get_count_fingerprint_mask (fps, remove_low_abundant = TRUE,
                                        remove_single_value = TRUE)
    expect_equal(mask, c(T, T, F, T))
})


context("Binary fingerprint mask")
test_that("is_correct", {
    fps <- matrix (c(
        F, F, T, F,
        F, T, T, T,
        F, F, T, F
    ), nrow = 3, byrow = TRUE)

    mask <- get_fingerprint_mask (fps, remove_low_variance = FALSE)
    expect_equal(mask, c(F, T, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_single_value = FALSE,
                                  remove_low_variance = FALSE)
    expect_equal(mask, c(T, T, T, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_redundant = FALSE)
    expect_equal(mask, c(F, T, F, T))
    mask <- get_fingerprint_mask (fps,
                                  remove_redundant = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, T, T, T))

    fps <- matrix (c(
        T, T, T, T,
        T, T, T, T,
        T, T, T, T
    ), nrow = 3, byrow = TRUE)

    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE)
    expect_equal(mask, c(F, F, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, F, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_redundant = FALSE)
    expect_equal(mask, c(F, F, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_redundant = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, T, T, T))


    mask <- get_fingerprint_mask (! fps,
                                  remove_low_variance = FALSE)
    expect_equal(mask, c(F, F, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, F, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_redundant = FALSE)
    expect_equal(mask, c(F, F, F, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_redundant = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, T, T, T))

    fps <- matrix (c(
        T, F, F, F,
        F, T, F, T,
        F, F, T, F
    ), nrow = 3, byrow = TRUE)

    mask <- get_fingerprint_mask (fps, remove_low_variance = FALSE)
    expect_equal(mask, c(T, T, T, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, T, T, F))
    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_redundant = FALSE)
    expect_equal(mask, c(T, T, T, T))

    fps <- matrix (c(
        T, F, F, T,
        F, T, F, T,
        F, F, T, F
    ), nrow = 3, byrow = TRUE)

    mask <- get_fingerprint_mask (fps, remove_low_variance = FALSE)
    expect_equal(mask, c(T, T, T, T))

    fps <- matrix (c(
        F, F, F, T, T,
        F, F, F, T, F,
        F, F, F, F, T
    ), nrow = 3, byrow = TRUE)

    mask <- get_fingerprint_mask (fps,
                                  remove_low_variance = FALSE,
                                  remove_single_value = FALSE)
    expect_equal(mask, c(T, F, F, T, T))
})

test_that ("border_cases_are_handled", {
    fps <- matrix(NA, nrow = 0, ncol = 0)
    mask <- get_fingerprint_mask(fps)
    expect_true(length(mask) == 0)

    fps <- matrix(c(T, T, F, T), nrow = 1)
    expect_equal(get_fingerprint_mask(fps), c(F, F, F, F))
    expect_equal(get_fingerprint_mask(fps, F, F, F), c(T, T, T, T))
    expect_equal(get_fingerprint_mask(fps, F, T, F), c(T, F, T, F))

    fps <- t(fps)
    expect_equal(get_fingerprint_mask(fps, T, T, F), T)
    expect_equal(get_fingerprint_mask(fps, T, T, T, 0.75), F)

    fps <- matrix(c(T, T, T, T), ncol = 1)
    expect_equal(get_fingerprint_mask(fps), F)
})

test_that ("wrong_input_is_handled", {
    fps <- matrix(1:3, ncol = 3)
    expect_error(get_fingerprint_mask(fps))
})

test_that ("low_variance_removal_is_correct", {
    fps <- matrix (c(
        T, F, F, T,
        F, T, T, T,
        F, T, T, F,
        F, F, T, F
    ), nrow = 4, byrow = TRUE)

    mask <- get_fingerprint_mask (fps, low_variance_tshd = 0.75,
                                  remove_single_value = FALSE, remove_redundant = FALSE)
    expect_equal(mask, c(F, T, F, T))

    mask <- get_fingerprint_mask (fps, low_variance_tshd = 0.5,
                                  remove_single_value = FALSE, remove_redundant = FALSE)
    expect_equal(mask, c(F, F, F, F))

    mask <- get_fingerprint_mask (fps, low_variance_tshd = 0.90,
                                  remove_single_value = FALSE, remove_redundant = FALSE)
    expect_equal(mask, c(T, T, T, T))
})

context ("Create set-difference fingerprints")
test_that ("Difference fingerprints are correct", {
    inchis <- c("InChI=1S/C9H10O4/c10-7-3-1-6(2-4-7)5-8(11)9(12)13/h1-4,8,10-11H,5H2,(H,12,13)",
                "InChI=1S/C4H8O3/c1-2-3(5)4(6)7/h3,5H,2H2,1H3,(H,6,7)")

    fps_4 <- calculate_fingerprints_from_inchi(
        inchis, fp_type = "circular", fp_mode = "count", circular.type="ECFP4")
    fps_6 <- calculate_fingerprints_from_inchi(
        inchis, fp_type = "circular", fp_mode = "count", circular.type="ECFP6")

    fps_diff <- setdiff_fingerprints(fps_6, fps_4)

    expect_equal(length(fps_4), length(fps_diff))
    expect_equal(length(fps_6), length(fps_diff))

    # Look at the individual fingerprint vectors
    fp_diff_1 <- fps_diff[[1]]
    fp_diff_2 <- fps_diff[[2]]

    fp_diff_feat_1 <- sapply(fp_diff_1@features, fingerprint::feature)
    fp_diff_feat_2 <- sapply(fp_diff_2@features, fingerprint::feature)

    expect_equal(fp_diff_feat_1, c("-1244114267", "-1028866103", "295780559", "722031324", "1602095716", "1972652850"))
    expect_equal(length(fp_diff_feat_2), 0)

    fp_diff_count_1 <- sapply(fp_diff_1@features, fingerprint::count)
    fp_diff_count_2 <- sapply(fp_diff_2@features, fingerprint::count)

    expect_equal(fp_diff_count_1, c(1, 1, 1, 1, 1, 1))
    expect_equal(length(fp_diff_count_2), 0)
})

context ("Create fingerprint matrix for hashed fps")
test_that ("Fingerprint matrix is correct", {
    make_mol_feat <- function(i) {
        new("feature", feature = LETTERS[i], count=as.integer(i))
    }

    fps <- list(
        "MOL1" = new("featvec", features = lapply(1:6, make_mol_feat)),
        "MOL2" = new("featvec", features = lapply(c(8,3,5), make_mol_feat)),
        "MOL3" = new("featvec", features = lapply(c(1,2,3,9,8), make_mol_feat)))

    fps_matrix <- fingerprints_to_matrix(fps, is_hashed = TRUE,
                                         sort_hash_keys = TRUE,
                                         add_colnames = TRUE)

    expect_equal(nrow(fps_matrix), length(fps))
    expect_equal(ncol(fps_matrix), 8)
    expect_equal(colnames(fps_matrix), LETTERS[c(1:6, 8, 9)])

    fps_matrix_ref <- matrix(c(1, 2, 3, 4, 5, 6, 0, 0,
                               0, 0, 3, 0, 5, 0, 8, 0,
                               1, 2, 3, 0, 0, 0, 8, 9),
                             nrow = length(fps), ncol = 8, byrow = TRUE)
    rownames(fps_matrix_ref) <- names(fps)
    colnames(fps_matrix_ref) <- LETTERS[c(1:6, 8, 9)]
    expect_equal(fps_matrix, fps_matrix_ref)

    # Do not sort the hash keys
    fps_matrix <- fingerprints_to_matrix(fps[c(2,3,1)], is_hashed = TRUE,
                                         sort_hash_keys = FALSE,
                                         add_colnames = TRUE)

    expect_equal(nrow(fps_matrix), length(fps))
    expect_equal(ncol(fps_matrix), 8)
    expect_equal(colnames(fps_matrix), LETTERS[c(8,3,5,1,2,9,4,6)])

    fps_matrix_ref <- matrix(c(8, 3, 5, 0, 0, 0, 0, 0,
                               8, 3, 0, 1, 2, 9, 0, 0,
                               0, 3, 5, 1, 2, 0, 4, 6),
                             nrow = length(fps), ncol = 8, byrow = TRUE)
    rownames(fps_matrix_ref) <- names(fps[c(2,3,1)])
    colnames(fps_matrix_ref) <- LETTERS[c(8,3,5,1,2,9,4,6)]
    expect_equal(fps_matrix, fps_matrix_ref)
})

context ("Fingerprints to JSON file")
test_that ("Fingerprints are put correctly converted to the JSON file", {
    fps <- list("MOL1" = new("featvec", features = lapply(1:6, function(i) {
            new("feature", feature = LETTERS[i], count=i)
         })),
         "MOL2" = new("featvec", features = lapply(7:12, function(i) {
             if(i %% 2 == 0){
                 c <- i
             } else {
                 c <- 0
             }
             new("feature", feature = LETTERS[i], count=as.integer(c))
         })))

    # $MOL1
    # Feature fingerprint
    # name =
    #     source =
    #     features =  A:1 B:2 C:3 D:4 E:5 F:6
    #
    # $MOL2
    # Feature fingerprint
    # name =
    #     source =
    #     features =  G:0 H:8 I:0 J:10 K:0 L:12

    tmp_file <- tempfile()
    write_fingerprint_to_json_file(fps, tmp_file)
    fps_json_string <- readLines(tmp_file)
    expect_equal(fps_json_string, '{"MOL1":{"A":[1],"B":[2],"C":[3],"D":[4],"E":[5],"F":[6]},"MOL2":{"G":[0],"H":[8],"I":[0],"J":[10],"K":[0],"L":[12]}}')

    tmp_file <- tempfile()
    write_fingerprint_to_json_file(fps, tmp_file, exclude_zero_fp = TRUE)
    fps_json_string <- readLines(tmp_file)
    expect_equal(fps_json_string, '{"MOL1":{"A":[1],"B":[2],"C":[3],"D":[4],"E":[5],"F":[6]},"MOL2":{"H":[8],"J":[10],"L":[12]}}')
})

test_that ("Binary fingerprints are correctly converted to the JSON file", {
    fps <- list(
        "MOL1" = new("fingerprint", nbit=8, bits=c(1,2,3,6,8)),
        "MOL2" = new("fingerprint", nbit=8, bits=c(1,7)))

    # $MOL1
    # Fingerprint object
    # name =
    #     length =  8
    # folded =  FALSE
    # source =
    #     bits on =  1 2 3 6 8
    #
    # $MOL2
    # Fingerprint object
    # name =
    #     length =  8
    # folded =  FALSE
    # source =
    #     bits on =  1 7

    tmp_file <- tempfile()
    write_fingerprint_to_json_file(fps, tmp_file)
    fps_json_string <- readLines(tmp_file)
    expect_equal(fps_json_string, '{"MOL1":{"1":[1],"2":[1],"3":[1],"4":[0],"5":[0],"6":[1],"7":[0],"8":[1]},"MOL2":{"1":[1],"2":[0],"3":[0],"4":[0],"5":[0],"6":[0],"7":[1],"8":[0]}}')

    tmp_file <- tempfile()
    write_fingerprint_to_json_file(fps, tmp_file, exclude_zero_fp = TRUE)
    fps_json_string <- readLines(tmp_file)
    expect_equal(fps_json_string, '{"MOL1":{"1":[1],"2":[1],"3":[1],"6":[1],"8":[1]},"MOL2":{"1":[1],"7":[1]}}')
})

test_that ("Unsupported fingerprint representation throws exception", {
    fps <- list("MOL1" = c(0,1,1,0,0), "MOL2" = c(0,1,1,1,1))
    expect_error(write_fingerprint_to_json_file(fps))
})

test_that ("Parameters are passed down to toJSON", {
    fps <- list("MOL1" = new("featvec", features = lapply(1:6, function(i) {
        new("feature", feature = LETTERS[i], count=i)
    })),
    "MOL2" = new("featvec", features = lapply(7:12, function(i) {
        if(i %% 2 == 0){
            c <- i
        } else {
            c <- 0
        }
        new("feature", feature = LETTERS[i], count=as.integer(c))
    })))

    tmp_file <- tempfile()
    write_fingerprint_to_json_file(fps, tmp_file, exclude_zero_fp = TRUE, auto_unbox = TRUE)
    fps_json_string <- readLines(tmp_file)
    expect_equal(fps_json_string, '{"MOL1":{"A":1,"B":2,"C":3,"D":4,"E":5,"F":6},"MOL2":{"H":8,"J":10,"L":12}}')

    tmp_file <- tempfile()
    write_fingerprint_to_json_file(fps, tmp_file, auto_unbox = TRUE, pretty = TRUE)
    fps_json_string <- readLines(tmp_file)
    expect_equal(paste0(fps_json_string, collapse = "\n"),
'{
  "MOL1": {
    "A": 1,
    "B": 2,
    "C": 3,
    "D": 4,
    "E": 5,
    "F": 6
  },
  "MOL2": {
    "G": 0,
    "H": 8,
    "I": 0,
    "J": 10,
    "K": 0,
    "L": 12
  }
}')
})

context ("Fingerprints to matrix")
test_that ("Counting fps are converted correctly", {
    fps <- list(
        "MOL0" = NULL,
        "MOL1" = new("featvec", features = lapply(1:6, function(i) {
            new("feature", feature = LETTERS[i], count=i)
        })),
        "MOL2" = new("featvec", features = lapply(7:12, function(i) {
            if(i %% 2 == 0){
                c <- i
            } else {
                c <- 0
            }
            new("feature", feature = LETTERS[i], count=as.integer(c))
        })))

    fps_mat <- fingerprints_to_matrix(fps)

    expect_equal(nrow(fps_mat), 3)
    expect_equal(ncol(fps_mat), 6)
    expect_equal(names(fps), rownames(fps_mat))

    expect_true(all(is.na(fps_mat[1,])))
    expect_equal(fps_mat[2,], 1:6)
    expect_equal(fps_mat[3,], c(0,8,0,10,0,12))
})

test_that ("Binary fps are converted correctly", {
    fps <- list(
        "MOL1" = new("fingerprint", nbit=8, bits=c(1,2,3,6,8)),
        "MOL2" = NULL,
        "MOLX" = NULL,
        "MOL3" = new("fingerprint", nbit=8, bits=c(1,7)))

    fps_mat <- fingerprints_to_matrix(fps)

    expect_equal(nrow(fps_mat), 4)
    expect_equal(ncol(fps_mat), 8)
    expect_equal(names(fps), rownames(fps_mat))

    expect_true(all(is.na(fps_mat[2,])))
    expect_true(all(is.na(fps_mat[3,])))
    expect_equal(fps_mat[1,], c(1,1,1,0,0,1,0,1))
    expect_equal(fps_mat[4,], c(1,0,0,0,0,0,1,0))
})

test_that ("Invalid input is handled correctly", {
    # all fingerprints are zero
    fps <- list("MOL1" = NULL, "MOL2" = NULL)
    expect_error(fingerprints_to_matrix(fps))

    # fingerprints have different dimension
    fps <- list(
        "MOL1" = new("fingerprint", nbit=8, bits=c(1,2,3,6,8)),
        "MOL2" = new("fingerprint", nbit=7, bits=c(1,7)))
    expect_error(fingerprints_to_matrix(fps))
})

context ("Write fingerprints to csv-file")
test_that ("Fps are correctly stored", {
    fps <- list(
        "MOL1" = new("fingerprint", nbit=8, bits=c(1,2,3,6,8)),
        "MOL2" = NULL,
        "MOLX" = NULL,
        "MOL3" = new("fingerprint", nbit=8, bits=c(1,7)))

    tmpfile = tempfile()
    write_fingerprint_to_csv_file(fps, tmpfile)

    fps_csv <- read.table(tmpfile, stringsAsFactors = FALSE, sep = ",", header = FALSE)

    expect_equal(nrow(fps_csv), length(fps))
    expect_equal(ncol(fps_csv), 1 + 8)
    expect_equal(fps_csv[,1], names(fps))
})

context ("Write fingerprint mask to file")
test_that("Mask is correctly stored", {
    mask <- matrix(c(T, T, F, F, T), nrow = 1)

    tmpfile = tempfile()
    write_fingerprint_mask_to_csv_file(mask, tmpfile)

    mask_csv <- read.table(tmpfile, sep = ",", header = FALSE)

    expect_true(all(as.matrix(mask_csv) == mask))
})

context ("Package data")
test_that ("MACCS_SMARTS_PATTERN_are_accessable", {
    expect_equal(length(count_maccs_pattern), 142)
})

context ("Molecule configuration specific for fingerprint type")
test_that ("MACCS_counting_fps_are_calculated", {
    fp = calculate_fingerprints_from_smiles("c1ccccc1", fp_type="maccs", fp_mode="count")
    expect_equal(length(fp[[1]]@features), length(count_maccs_pattern))
    expect_equal(fingerprint::count(fp[[1]]@features[[112]]), 1)
})
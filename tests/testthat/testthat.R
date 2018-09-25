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

    expect_true(! is.null(fps[[1]]))
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

    expect_true(! is.null(fps[[1]]))
    expect_true(is.null(fps[[2]]))
    expect_true(! is.null(fps[[3]]))
    expect_true(is.null(fps[[4]]))
    expect_true(! is.null(fps[[5]]))
    expect_true(! is.null(fps[[6]]))
    expect_true(! is.null(fps[[7]]))
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

context("Descriptor matrix")
test_that("Descriptor matrix is correctly created", {
    smiles <- c("CCC(C)Cl", "NC(C)C(=Oasd)O", "F/C=C/Fasdasd", "CCC(C)Cl")

    # Check descriptor with a single value
    desc <- calculate_descriptor_from_smiles(
        smiles, desc_name=rcdk::get.desc.names()[4])

    desc_mat <- descriptors_to_matrix(desc)

    expect_equal(rownames(desc_mat), smiles)
    expect_equal(nrow(desc_mat), 4)
    expect_equal(ncol(desc_mat), 1)
    expect_false(is.na(desc_mat[1,]))
    expect_true(is.na(desc_mat[2,]))
    expect_true(is.na(desc_mat[3,]))
    expect_false(is.na(desc_mat[4,]))

    # Check descriptor with a single value
    desc <- calculate_descriptor_from_smiles(
        smiles, desc_name=rcdk::get.desc.names()[2])

    desc_mat <- descriptors_to_matrix(desc)

    expect_equal(rownames(desc_mat), smiles)
    expect_equal(nrow(desc_mat), 4)
    expect_equal(ncol(desc_mat), 11)
    expect_true(! any(is.na(desc_mat[1,])))
    expect_true(all(is.na(desc_mat[2,])))
    expect_true(all(is.na(desc_mat[3,])))
    expect_true(! any(is.na(desc_mat[4,])))
})
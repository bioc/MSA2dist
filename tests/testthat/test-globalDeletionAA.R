test_that("globalDeletionAA()", {
    cds1 <- Biostrings::DNAString("ATGCAACATTGC")
    cds2 <- Biostrings::DNAString("ATG---CATTGC")
    cds1.cds2.aln <- c(Biostrings::DNAStringSet(cds1),
        Biostrings::DNAStringSet(cds2))
    aa.aln <- cds1.cds2.aln |> cds2aa()
    expect_true(as.character(globalDeletionAA(aa.aln)[1]) == "MHC")
})

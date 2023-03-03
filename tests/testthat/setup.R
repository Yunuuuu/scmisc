set.seed(1L)
sce <- scuttle::mockSCE(ncells = 800L)
sce$batch <- paste0(
    "sample",
    sample(1:4, size = 800L, replace = TRUE)
)
colnames(sce) <- paste(sce$batch, colnames(sce), sep = "_")
sce$block <- paste0("block", sample(1:2, size = 800L, replace = TRUE))

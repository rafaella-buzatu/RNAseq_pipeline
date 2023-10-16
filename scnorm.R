library(SCnorm)
pathToData = 'data/Sample_MUC30433/filteredCountMatrix.csv'

countMatrix = as.matrix(read.csv(pathToData, row.names = 1))

scData <- SingleCellExperiment::SingleCellExperiment(assays = list('counts' = countMatrix) )


DataNorm <- SCnorm(Data = scData,  PrintProgressPlots = TRUE, Conditions = rep("1", ncol(scData)), FilterCellNum = 10, K = 1, NCores=40, reportSF = TRUE)

NormalizedData <- SingleCellExperiment::normcounts(DataNorm)

library(purrr)
path_to_output <- "data/"
path_to_datasets <- "data/chunks/"
datasets <- dir(path_to_datasets)
chunks_by_dataset <- map(datasets, ~ paste0(.x,"/",dir(paste0(path_to_datasets,.x))))
names(chunks_by_dataset) <- datasets

KFold <- function(s,k) {
    len_s <- length(s)
    if (k > len_s) {
        k <- len_s
    }
    s <- s[sample(len_s)]
    step <- as.integer(len_s/k)
    remainder <- len_s%%k
    folds <- map(seq(1,len_s - remainder,step), ~ s[.x:(.x+step-1)])
    if (remainder!=0) {
        for (i in 1:remainder) 
        {
            folds[[i]] <- c(folds[[i]], s[[len_s-i+1]])
        }
    }
    folds
}

dp.peaks.sets <- modify(chunks_by_dataset, ~ KFold(.x,10))
save(dp.peaks.sets, file=paste0(path_to_output,"dp.peaks.sets.RData"))
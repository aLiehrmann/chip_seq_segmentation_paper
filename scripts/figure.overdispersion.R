library(data.table)
library(ggplot2)
library(parallel)
library(purrr)


path.files <- "data/chunks/"

load("data/res_cv_PDPA_gaussian_linear_lambda_10.RData")
load("data/all_modelSelection_PDPA.RData")

l.gaussian <- list()
z = 1
t.files <- unique(all.modelSelection[,list(chunk.name)])
index.current.file <- 1
l.gaussian  <- mclapply(1:nrow(t.files), function(i.file)
{
    print(paste0(i.file,"/",nrow(t.files)))
    index.current.file <- index.current.file + 1
    current.file <- t.files[i.file,]$chunk.name
    load(paste0(path.files,current.file,"/counts.RData"))
    load(paste0(path.files,current.file,"/PDPA.RData"))
    t.samples <- res_models[chunk.name == current.file]
    return(do.call(rbind, lapply(1:nrow(t.samples), function(i.sample)
    {
        current.sample <- t.samples[i.sample,]$sample.id
        df.current.counts <- counts[
            counts$sample.id==current.sample,
        ]

        cp_star <- res_models[
            sample.id == current.sample &
            chunk.name == current.file
        ]$changepoints

        t.current.model <- dt.models[
            sample == current.sample
            & changepoints == cp_star
        ]

        data <- rep(
            df.current.counts$coverage, 
            times = df.current.counts$chromEnd - df.current.counts$chromStart
        )

        cp <- c(1,t.current.model$fit[[1]]$changepoints)
        
        mu <- unlist(lapply(
            2:length(cp), 
            function(x){
                mean(rep(
                    sqrt(df.current.counts$coverage[cp[x-1]:cp[x]]+3/8),
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        var <- unlist(lapply(
            2:length(cp), 
            function(x){
                var(rep(
                    sqrt(df.current.counts$coverage[cp[x-1]:cp[x]]+3/8),
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        mu.raw <- unlist(lapply(
            2:length(cp), 
            function(x){
                mean(rep(
                    df.current.counts$coverage[cp[x-1]:cp[x]],
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        var.raw <- unlist(lapply(
            2:length(cp), 
            function(x){
                var(rep(
                    df.current.counts$coverage[cp[x-1]:cp[x]],
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        cp <- c(0, t.current.model$fit[[1]]$changepoints)
        weights <- df.current.counts$chromEnd - df.current.counts$chromStart
        counts.raw <- df.current.counts$coverage
        counts <- sqrt(df.current.counts$coverage+3/8)
        vec_mu <- rep(mu, times = diff(cp))
        vec_mu.raw <- rep(mu.raw, times = diff(cp))
        v <- sum((counts - vec_mu)^2 * weights) / (sum(weights) -(length(cp)-1))
        v.raw <- sum((counts.raw - vec_mu.raw)^2 * weights) / (sum(weights)-(length(cp)-1))
        
        z <- z+1

        return(data.table(
            id = z,
            empirical_mean = mu.raw,
            theoretical_variance = v,
            empirical_variance = var,
            type = strsplit(current.file,"_")[[1]][1]
        ))


    })))
}, mc.cores = 8)



l.gaussian  <- do.call(rbind, l.gaussian )

l.gaussian [,
    type:=ifelse(
        type == "H3K36me3",
        "H3K36me3\n(broad peaks)",
        "H3K4me3\n(sharp peaks)"
    )
]

l.gaussian [,model:="gaussian"]

load("data/res_cv_updown_poisson_linear_lambda_4.MODELS.RData")
load("data/all_modelSelection_updown_poisson.RData")


l.poisson <- list()
z = 1
t.files <- unique(all.modelSelection[,list(chunk.name)])
index.current.file <- 1
l.poisson <- mclapply(1:nrow(t.files), function(i.file)
{
    print(paste0(i.file,"/",nrow(t.files)))
    index.current.file <- index.current.file + 1
    current.file <- t.files[i.file,]$chunk.name
    load(paste0(path.files,current.file,"/counts.RData"))
    load(paste0(path.files,current.file,"/updown_poisson.RData"))
    t.samples <- res_models[chunk.name == current.file]
    return(do.call(rbind, lapply(1:nrow(t.samples), function(i.sample)
    {
        current.sample <- t.samples[i.sample,]$sample.id
        df.current.counts <- counts[
            counts$sample.id==current.sample,
        ]

        cp_star <- res_models[
            sample.id == current.sample &
            chunk.name == current.file
        ]$changepoints

        t.current.model <- dt.models[
            sample == current.sample
            & changepoints == cp_star
        ]

        data <- rep(
            df.current.counts$coverage, 
            times = df.current.counts$chromEnd - df.current.counts$chromStart
        )

        cp <- c(1,t.current.model$fit[[1]]$changepoints)
        mu <- unlist(lapply(
            2:length(cp), 
            function(x){
                mean(rep(
                    df.current.counts$coverage[cp[x-1]:cp[x]],
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        var <- unlist(lapply(
            2:length(cp), 
            function(x){
                var(rep(
                    df.current.counts$coverage[cp[x-1]:cp[x]],
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        z <- z+1

        return(data.table(
            id = z,
            empirical_mean = mu,
            theoretical_variance = mu,
            empirical_variance = var,
            type = strsplit(current.file,"_")[[1]][1]
        ))


    })))
}, mc.cores = 8)



l.poisson <- do.call(rbind, l.poisson)

l.poisson[,
    type:=ifelse(
        type == "H3K36me3",
        "H3K36me3\n(broad peaks)",
        "H3K4me3\n(sharp peaks)"
    )
]

l.poisson[,model:="poisson"]


load("data/res_cv_PDPA_negbin_linear_lambda_constant_phi_4.MODELS.RData")
load("data/all_modelSelection_PDPA_negbin_2.RData")

l.negbin <- list()
z = 1
t.files <- unique(all.modelSelection[,list(chunk.name)])
index.current.file <- 1
l.negbin <- mclapply(1:nrow(t.files), function(i.file)
{
    print(paste0(i.file,"/",nrow(t.files)))
    index.current.file <- index.current.file + 1
    current.file <- t.files[i.file,]$chunk.name
    load(paste0(path.files,current.file,"/counts.RData"))
    load(paste0(path.files,current.file,"/PDPA_negbin_2.RData"))
    t.samples <- res_models[chunk.name == current.file]
    return(do.call(rbind, lapply(1:nrow(t.samples), function(i.sample)
    {
        current.sample <- t.samples[i.sample,]$sample.id
        df.current.counts <- counts[
            counts$sample.id==current.sample,
        ]

        cp_star <- res_models[
            sample.id == current.sample &
            chunk.name == current.file
        ]$changepoints

        phi_star <- res_models[
            sample.id == current.sample &
            chunk.name == current.file
        ]$phi

        t.current.model <- dt.models[
            sample == current.sample
            & changepoints == cp_star
            & phi == phi_star
        ]

        data <- rep(
            df.current.counts$coverage, 
            times = df.current.counts$chromEnd - df.current.counts$chromStart
        )

        cp <- c(1,t.current.model$fit[[1]]$changepoints)
        mu <- unlist(lapply(
            2:length(cp), 
            function(x){
                mean(rep(
                    df.current.counts$coverage[cp[x-1]:cp[x]],
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        var <- unlist(lapply(
            2:length(cp), 
            function(x){
                var(rep(
                    df.current.counts$coverage[cp[x-1]:cp[x]],
                    times = df.current.counts$chromEnd[cp[x-1]:cp[x]] - df.current.counts$chromStart[cp[x-1]:cp[x]]
                ))
            }
        ))

        z <- z+1

        return(data.table(
            id = z,
            empirical_mean = mu,
            empirical_variance = var,
            phi_learned = phi_star,
            type = strsplit(current.file,"_")[[1]][1]
        ))


    })))
}, mc.cores = 8)



l.negbin <- do.call(rbind, l.negbin)

l.negbin[,
    type:=ifelse(
        type == "H3K36me3",
        "H3K36me3\n(broad peaks)",
        "H3K4me3\n(sharp peaks)"
    )
]

l.negbin[,theoretical_variance:=empirical_mean+(1/(phi_learned))*empirical_mean^2]
l.negbin[,model:="negative binomial"]
l.negbin[,phi_learned:=NULL]

l <- rbind(l.poisson,l.negbin, l.gaussian)



l$model <- factor(l$model, levels = c("poisson","negative binomial", "gaussian"))


lab.mean <- data.table(
    label = c(
        "segments with < 1\ncounts in average",
        "segments with > 1\ncounts in average",
        "segments with < 1\ncounts in average",
        "segments with > 1\ncounts in average"
    ),
    x = c(-3.2, 3.2, -3.5, 3.5),
    y = c(-11),
    type = c(
       "H3K36me3\n(broad peaks)",
       "H3K36me3\n(broad peaks)",
       "H3K4me3\n(sharp peaks)",
       "H3K4me3\n(sharp peaks)"
    )
)

lab.model <- data.table(
    label = "model assumption",
    y = c(1.5),
    x = c(-10.6, -7.75),
    type = c(
       "H3K36me3\n(broad peaks)",
       "H3K4me3\n(sharp peaks)"
    )
)

l[,x := log2(empirical_mean)]
l[,y := log2(empirical_variance/theoretical_variance)]


median_models <- l[, list(y=median(y, na.rm = TRUE)),by=list(type, model)]
l[model=="negative binomial"]$model <- "negative\nbinomial"
l$model <- factor(l$model, , levels = c("poisson","negative\nbinomial", "gaussian"))

g6 <- ggplot(
    aes(
    x=factor(model), 
    y=y),
    data = l
)+
facet_grid(.~type)+
geom_hline(yintercept=0, color="red")+
geom_violin(alpha=0.1)+
geom_boxplot(width=0.2, fill=0.1)+
ylab("log2(empirical variance / theoretical variance)")+
xlab("model")+
labs(color= "model")+
theme_bw()+
theme(
    legend.position="bottom",
    text=element_text(size=18),
    strip.background = element_rect(fill="grey95"),
    axis.text.x=element_text(angle=45,hjust=1),
    panel.spacing = unit(0, "lines")
)

g7 <- g6+coord_cartesian(ylim=c(-2.5,3))+ylab(NULL)

library("gridExtra")
grid.arrange(g6, g7, ncol=2)

pdf("figures/figure_modeling_of_overdispersion.pdf",
  width = 12,
  height = 6
)
grid.arrange(g6, g7, ncol=2)
dev.off()
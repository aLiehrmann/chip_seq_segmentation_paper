library(data.table)
library(ggplot2)
library(furrr)
plan(multiprocess, workers=10)


load("data/all_modelSelection_PDPA_negbin.RData")
load("data/dp.peaks.sets.RData")
setkey(all.modelSelection,chunk.name)
training_chunks <- unlist(dp.peaks.sets$H3K4me3_TDH_other[-4])
validation_chunks <- dp.peaks.sets$H3K4me3_TDH_other[[4]]
training_models <- all.modelSelection[training_chunks]
validation_models <- all.modelSelection[validation_chunks]
training_models <- training_models[!is.na(min.log.lambda),]
validation_models <- validation_models[!is.na(min.log.lambda),]
training_models <- split(training_models, by="phi")
dt <- do.call(rbind, future_map(training_models, function(dt_by_phi) {
    sorted.lambda <- sort(c(dt_by_phi$min.log.lambda, dt_by_phi$max.log.lambda))
    grid <- c(sorted.lambda[sorted.lambda[2:length(sorted.lambda)]!=sorted.lambda[1:(length(sorted.lambda)-1)]],sorted.lambda[length(sorted.lambda)])
    errors <- vector(mode="integer", length= length(grid)-1)
    walk(1:nrow(dt_by_phi), function(i){
      a <- which(grid>=dt_by_phi$min.log.lambda[[i]] & grid<dt_by_phi$max.log.lambda[[i]])
      errors[a] <<- errors[a] + dt_by_phi$total.errors[[i]]
    })
    data.table(min.log.lambda = grid[1:(length(grid)-1)], max.log.lambda = grid[2:(length(grid))], errors = errors, log.phi = log(dt_by_phi$phi[[1]]))
}, .options = furrr_options(seed=2020)))


predicted.phi <- dt[which.min(errors),]$log.phi
tmp <- dt[errors==min(errors),]
predicted.lambda <- (tmp[1,]$min.log.lambda + tmp[nrow(tmp),]$max.log.lambda) / 2



minimum_by_phi <- do.call(rbind, lapply(split(dt,by="log.phi"), function(x){
    data.table(
        log.phi = x$log.phi[[1]],
        errors = min(x$errors)
    )
}))


library(colorspace)
library(latex2exp)

label_minimum_error <- data.table(

)

g1 <- ggplot(
    data = dt,
    aes(xmin=log.phi- 0.2, ymin=min.log.lambda, xmax=log.phi + 0.2, ymax=max.log.lambda)
    ) +
    geom_rect(aes(fill=errors, color=errors)) +
    geom_rect(
        xmin = predicted.phi - 0.2, 
        xmax = predicted.phi + 0.2, 
        ymin = predicted.lambda - 0.1, 
        ymax = predicted.lambda + 0.1, 
        fill="red"
    )+
    annotate(
        geom = "text",
        x = predicted.phi+0.1,
        y = predicted.lambda - 1.5,
        label = TeX("$\\log(\\phi^*)$"),
        color = "red",
        size= 6
    )+
    annotate(
        geom = "text",
        x = predicted.phi- 1.1,
        y = predicted.lambda,
        label = TeX("$\\log(\\lambda^*)$"),
        color = "red",
        size= 6
    )+
    scale_fill_continuous_sequential(palette = "Blues")+
    scale_color_continuous_sequential(palette = "Blues")+
    xlab(TeX("dispersion: $\\log(\\phi)$"))+
    ylab(TeX("penalty: $\\log(\\lambda)$"))+
    theme_bw()+
    theme(
        text = element_text(size=20),
    )+
    labs(fill="errors:")+
    scale_x_continuous(breaks=round(seq(log(1),log(10000),length=16),1))+
    guides(color=FALSE)


g2 <- ggplot(aes(x=log.phi, y=errors), data=minimum_by_phi)+
    geom_line()+
    theme_bw()+
    geom_point(x=predicted.phi, y=min(minimum_by_phi$error), color="red", size=3)+
    theme(
        text = element_text(size=20),
    ) +
    xlab(TeX("dispersion: $\\log(\\phi)$"))+
    ylab(TeX("minimum errors"))+
    scale_x_continuous(breaks=round(seq(log(1),log(10000),length=16),1))

library("cowplot")
g <- ggdraw() +
    draw_plot(g1, x = 0, y = .5, width = 1, height = .5)+
    draw_plot(g2, x = 0, y = 0, width = 0.87, height = .5)

pdf(file="figures/figure_target_phi_constant_lambda_constant.pdf",
         height = 7,
          width=9)    
g
dev.off()
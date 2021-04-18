library(data.table)
library(purrr)
library(ggplot2)
library(stringr)

res_cv <- rbindlist(map(Sys.glob("data/test_error_vs_nb_labels*"), function(f){
  load(f)
  res_cv
}))

res_cv[, model := str_replace(res_cv$model, "PDPA", "unconstrained")]
res_cv[, typePeak := str_extract(pattern=".*?(?=_)",set.name)]
res_cv[, model := str_extract(pattern=".*?(?=_constant.*)", model)]
g_accuracy <- ggplot(
    data=res_cv[samples %in% c(1, 5, 10, 20, 30, 40),], 
    aes(x=factor(samples), y=accuracy)
)+
facet_grid(typePeak~model,scales="free")+
geom_boxplot()+
ylab("test accuracy %")+
xlab("number of genomic windows used for training")+
theme_bw()+
theme(
  strip.background = element_rect(fill="grey95"), 
  panel.margin=grid::unit(0, "lines"),
  text = element_text(size=13),
  legend.position = "bottom"
)

pdf(file = "figures/figure_test_error_vs_labels.pdf", width=7, height=15)
g_accuracy
dev.off()

pdf(file = "figures/figure_test_error_vs_labels.pdf", width=16, height=7)
g_accuracy
dev.off()


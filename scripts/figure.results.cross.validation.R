library(stringr)
library(data.table)
library(ggplot2)
library(tibble)
library(purrr)

signi <- function(p.values) {
  unlist(lapply(p.values, function(p.val){
    if (p.val > 0.05) {
      "ns"
    } else if (p.val > 0.01) {
      "*"
    } else if (p.val > 0.001) {
      "**"
    } else {
      "***"
    }
  }))
}

load("data/res_cv_all.RData")


tmp_macs <- res_cv_all[grep("H3K4me3", set.name)][model=="MACS"]
tmp_macs[,model := "HMCan (H3K36me3)\nMACS (H3K4me3)"]
tmp_hmcan <- res_cv_all[grep("H3K36me3", set.name)][model=="HMCan"]
tmp_hmcan[,model := "HMCan (H3K36me3)\nMACS (H3K4me3)"]
res_cv_all <- rbind(res_cv_all[model != "MACS" & model != "HMCan"], tmp_hmcan, tmp_macs)
res_cv_all[,type.peak := str_extract(res_cv_all$set.name, ".+3")]


###############################
## linear VS constant lambda ##
###############################


res_cv_all_linear_vs_constant_lambda <- res_cv_all[model != "HMCan (H3K36me3)\nMACS (H3K4me3)"]
res_cv_all_linear_vs_constant_lambda <- res_cv_all_linear_vs_constant_lambda[grep(res_cv_all_linear_vs_constant_lambda$model, pattern="linear_phi",invert=T),]

model.lab_ <- unlist(lapply(str_split(pattern="_", res_cv_all_linear_vs_constant_lambda$model), function(x){
  noise <- ifelse(x[[2]] == "negbin", "negative binomial", x[[2]])
  if (x[[1]] == "updown"){
    paste0(x[[1]]," ",noise)
  }
  else
  {
    rule <- ifelse(x[[3]] == "smallest", "thinnest", ifelse(x[[3]] == "largest", "largest","max jump"))
    paste0("unconstrained"," ",noise," ",rule)#paste0(x[[1]]," ",noise," ",rule)
  }
}))

model.lab2_ <- str_extract(pattern="(constant|linear)_lambda", res_cv_all_linear_vs_constant_lambda$model)
res_cv_all_linear_vs_constant_lambda[,model.lab := model.lab_]
res_cv_all_linear_vs_constant_lambda[,model.lab2 := model.lab2_]

dt_best <- do.call(rbind, lapply(split(res_cv_all_linear_vs_constant_lambda, by=c("model.lab","set.name")), function(x){
  res <- x[,mean(accuracy), by=list(model.lab2)]
  m <- res[which.max(V1),]$model.lab2
  data.table(
    model.lab = x$model.lab[[1]], 
    set.name = x$set.name[[1]],
    Best = paste0(str_split(pattern="_",m)[[1]][[1]]," lambda"),
    accuracy = res[which.max(V1),]$V1,
    model.lab2 = m
    )
}))

g <- ggplot(aes(y=as.factor(model.lab2), x=accuracy), data=res_cv_all_linear_vs_constant_lambda)+
  theme_bw()+
  stat_summary(fun.data = mean_se, geom="point", color="red")+
  stat_summary(fun.data = mean_se, geom="errorbar",color="red",width=0.5, size=0.5)+
  facet_grid(model.lab~set.name, scales ="free")+
  geom_rect(aes(xmin=-Inf,xmax=+Inf,ymin=-Inf, ymax=+Inf, fill=Best), data = dt_best, alpha=0.2)+
  #geom_rect(aes(xmin=-Inf,xmax=+Inf,ymin=-Inf, ymax=+Inf, fill=diff),alpha=0.2, data = lab_diff)+
  theme(
    strip.text.y = element_text(angle = 0),
    strip.background = element_rect(fill="grey95"),
    panel.margin=grid::unit(0, "lines")
  )+
  scale_fill_manual(values=c("#E69F00","#56B4E9"))+
  xlab("accuracy % on test folds")+
  ylab("models")+
  labs(fill = "Better in average")


pdf(file = "figures/figure_comparison_learning_methods.pdf",
  width = 17,
  height = 8,
)
print(g)
dev.off()

signi_linear_vs_constant_lambda <- do.call(rbind, map(split(res_cv_all_linear_vs_constant_lambda, by=c("model.lab","type.peak")), function(x){
  label <- mean_cl_normal(x[model.lab2 == "constant_lambda"]$accuracy - x[model.lab2 == "linear_lambda"]$accuracy)
  res <- t.test(x[model.lab2 == "constant_lambda"]$accuracy, x[model.lab2 == "linear_lambda"]$accuracy, paired=T)
  p.value <- res$p.value
  data.table(
    diff.accuracy = NaN,
    type.peak = x$type.peak[[1]],
    lab = x$model.lab[[1]],
    p.value = p.value,
    x = label$y,
    xmin = label$ymin,
    xmax = label$ymax
  )
}))


signi_linear_vs_constant_lambda [, adj.p.value := p.adjust(p.value, method="fdr")]
signi_linear_vs_constant_lambda [,label := signi(adj.p.value)]


#######################
## max jump vs other ##
#######################


# max_jump_vs_other <- tribble(
#   ~target,                                                          ~ref,                                           ~type.peak,   ~lab,                       ~lab2,
#   "PDPA_gaussian_largest_peak_rule_constant_lambda",                "PDPA_gaussian_constant_lambda",                "H3K36me3",   "PDPA gaussian",            "largest peak",
#   "PDPA_gaussian_smallest_peak_rule_constant_lambda",               "PDPA_gaussian_constant_lambda",                "H3K36me3",   "PDPA gaussian",            "thinnest peak",
#   "PDPA_poisson_largest_peak_rule_constant_lambda",                 "PDPA_poisson_constant_lambda",                 "H3K36me3",   "PDPA poisson",             "largest peak",
#   "PDPA_poisson_smallest_peak_rule_constant_lambda",                "PDPA_poisson_constant_lambda",                 "H3K36me3",   "PDPA poisson",             "thinnest peak",
#   "PDPA_negbin_largest_peak_rule_constant_lambda_constant_phi",     "PDPA_negbin_constant_lambda_constant_phi",     "H3K36me3",   "PDPA negative binomial",   "largest peak",
#   "PDPA_negbin_smallest_peak_rule_constant_lambda_constant_phi",    "PDPA_negbin_constant_lambda_constant_phi",     "H3K36me3",   "PDPA negative binomial",   "thinnest peak",
#   "PDPA_gaussian_largest_peak_rule_linear_lambda",                  "PDPA_gaussian_linear_lambda",                  "H3K4me3",    "PDPA gaussian",            "largest peak",
#   "PDPA_gaussian_smallest_peak_rule_linear_lambda",                 "PDPA_gaussian_linear_lambda",                  "H3K4me3",    "PDPA gaussian",            "thinnest peak",
#   "PDPA_poisson_largest_peak_rule_linear_lambda",                   "PDPA_poisson_linear_lambda",                   "H3K4me3",    "PDPA poisson",             "largest peak",
#   "PDPA_poisson_smallest_peak_rule_linear_lambda",                  "PDPA_poisson_linear_lambda",                   "H3K4me3",    "PDPA poisson",             "thinnest peak",
#   "PDPA_negbin_largest_peak_rule_linear_lambda_constant_phi",       "PDPA_negbin_linear_lambda_constant_phi",       "H3K4me3",    "PDPA negative binomial",   "largest peak",
#   "PDPA_negbin_smallest_peak_rule_linear_lambda_constant_phi",      "PDPA_negbin_linear_lambda_constant_phi",       "H3K4me3",    "PDPA negative binomial",   "thinnest peak"
# )

max_jump_vs_other <- tribble(
  ~target,                                                          ~ref,                                           ~type.peak,   ~lab,                       ~lab2,
  "PDPA_gaussian_largest_peak_rule_constant_lambda",                "PDPA_gaussian_constant_lambda",                "H3K36me3",   "unconstrained\ngaussian",            "largest peak",
  "PDPA_gaussian_smallest_peak_rule_constant_lambda",               "PDPA_gaussian_constant_lambda",                "H3K36me3",   "unconstrained\ngaussian",            "thinnest peak",
  "PDPA_poisson_largest_peak_rule_constant_lambda",                 "PDPA_poisson_constant_lambda",                 "H3K36me3",   "unconstrained\npoisson",             "largest peak",
  "PDPA_poisson_smallest_peak_rule_constant_lambda",                "PDPA_poisson_constant_lambda",                 "H3K36me3",   "unconstrained\npoisson",             "thinnest peak",
  "PDPA_negbin_largest_peak_rule_constant_lambda_constant_phi",     "PDPA_negbin_constant_lambda_constant_phi",     "H3K36me3",   "unconstrained\nnegative binomial",   "largest peak",
  "PDPA_negbin_smallest_peak_rule_constant_lambda_constant_phi",    "PDPA_negbin_constant_lambda_constant_phi",     "H3K36me3",   "unconstrained\nnegative binomial",   "thinnest peak",
  "PDPA_gaussian_largest_peak_rule_linear_lambda",                  "PDPA_gaussian_linear_lambda",                  "H3K4me3",    "unconstrained\ngaussian",            "largest peak",
  "PDPA_gaussian_smallest_peak_rule_linear_lambda",                 "PDPA_gaussian_linear_lambda",                  "H3K4me3",    "unconstrained\ngaussian",            "thinnest peak",
  "PDPA_poisson_largest_peak_rule_linear_lambda",                   "PDPA_poisson_linear_lambda",                   "H3K4me3",    "unconstrained\npoisson",             "largest peak",
  "PDPA_poisson_smallest_peak_rule_linear_lambda",                  "PDPA_poisson_linear_lambda",                   "H3K4me3",    "unconstrained\npoisson",             "thinnest peak",
  "PDPA_negbin_largest_peak_rule_linear_lambda_constant_phi",       "PDPA_negbin_linear_lambda_constant_phi",       "H3K4me3",    "unconstrained\nnegative binomial",   "largest peak",
  "PDPA_negbin_smallest_peak_rule_linear_lambda_constant_phi",      "PDPA_negbin_linear_lambda_constant_phi",       "H3K4me3",    "unconstrained\nnegative binomial",   "thinnest peak"
)

res_cv_all_max_jump_vs_other <- do.call(rbind,pmap(max_jump_vs_other, function(target, ref, type.peak, lab, lab2){
  type.peak_ <- type.peak
  cbind(
    res_cv_all[model == target & type.peak == type.peak_], 
    diff.accuracy = res_cv_all[model == target & type.peak == type.peak_]$accuracy - res_cv_all[model == ref & type.peak == type.peak_]$accuracy,
    lab = lab,
    lab2 = lab2
  )
}))

# target = "PDPA_gaussian_largest_peak_rule_constant_lambda"
# ref = "PDPA_gaussian_constant_lambda"
# type.peak = "H3K36me3"
# lab = "PDPA gaussian"
# lab2 = "largest peak"

target = "PDPA_gaussian_largest_peak_rule_constant_lambda"
ref = "PDPA_gaussian_constant_lambda"
type.peak = "H3K36me3"
lab = "unconstrained\ngaussian"
lab2 = "largest peak"

signi_max_jump_vs_other <- do.call(rbind,pmap(max_jump_vs_other, function(target, ref, type.peak, lab, lab2){
  type.peak_ <- type.peak
  label <- mean_cl_normal(res_cv_all[model == target & type.peak==type.peak_]$accuracy - res_cv_all[model == ref & type.peak==type.peak_]$accuracy)
  res <- t.test(res_cv_all[model == target & type.peak==type.peak_]$accuracy, res_cv_all[model == ref & type.peak==type.peak_]$accuracy, paired=T)
  p.value <- res$p.value
  data.table(
    diff.accuracy = NaN,
    type.peak = type.peak,
    lab = lab,
    lab2 = lab2,
    p.value = p.value,
    x = label$y,
    xmin = label$ymin,
    xmax = label$ymax
  )
}))

signi_max_jump_vs_other[, adj.p.value := p.adjust(p.value, method="fdr")]
signi_max_jump_vs_other[,label := signi(adj.p.value)]


signi_max_jump_vs_other[,type.peak := ifelse(type.peak=="H3K36me3","H3K36me3\n(broad peaks)","H3K4me3\n(sharp peaks)" )]
res_cv_all_max_jump_vs_other[,type.peak := ifelse(type.peak=="H3K36me3","H3K36me3\n(broad peaks)","H3K4me3\n(sharp peaks)" )]

g <- ggplot(aes(y=as.factor(lab2), x=diff.accuracy), data=res_cv_all_max_jump_vs_other)+
  theme_bw()+
  facet_grid(lab~type.peak, scales ="free")+
  geom_vline(xintercept = 0, color="Royalblue")+
  geom_errorbar(aes(xmin = xmin, xmax=xmax), data= signi_max_jump_vs_other, color="red",width=0.5, size=1)+
  geom_point(aes(x = x), data= signi_max_jump_vs_other, color="red")+
  geom_text(aes(x = xmin -1, y=as.factor(lab2), label=label), data= signi_max_jump_vs_other)+
  xlab("accuracy(target model)% - accuracy(reference model)% on test folds")+
  ylab("target model")+
  theme(
    strip.background = element_rect(fill="grey95"), 
    strip.text.y = element_text(angle = 0),
    panel.margin=grid::unit(0, "lines"),
    text = element_text(size=15)
  )+
  ggtitle("reference model : max jump post-processing rule")

pdf(file = "figures/figure_comparison_post_processing_tules.pdf",
  width = 8,
  height = 3,
)
print(g)
dev.off()


#####################
## updown vs other ##
#####################

# updown_poisson_vs_other <- tribble(
#   ~target,                                        ~ref,                              ~type.peak,   ~lab2,                      ~lab,
#   "PDPA_gaussian_constant_lambda",                "updown_poisson_constant_lambda",  "H3K36me3",   "gaussian",                 "PDPA",
#   "PDPA_poisson_constant_lambda",                 "updown_poisson_constant_lambda",  "H3K36me3",   "poisson",                  "PDPA",
#   "PDPA_negbin_constant_lambda_constant_phi",     "updown_poisson_constant_lambda",  "H3K36me3",   "negative binomial",        "PDPA",
#   "updown_gaussian_constant_lambda",              "updown_poisson_constant_lambda",  "H3K36me3",   "gaussian",                 "updown",
#   "updown_negbin_constant_lambda_constant_phi",   "updown_poisson_constant_lambda",  "H3K36me3",   "negative binomial",        "updown",
#   "PDPA_gaussian_linear_lambda",                  "updown_poisson_linear_lambda",    "H3K4me3",    "gaussian",                 "PDPA",
#   "PDPA_poisson_linear_lambda",                   "updown_poisson_linear_lambda",    "H3K4me3",    "poisson",                  "PDPA",
#   "PDPA_negbin_linear_lambda_constant_phi",       "updown_poisson_linear_lambda",    "H3K4me3",    "negative binomial",        "PDPA",
#   "updown_gaussian_linear_lambda",                "updown_poisson_linear_lambda",    "H3K4me3",    "gaussian",                 "updown",
#   "updown_negbin_linear_lambda_constant_phi",     "updown_poisson_linear_lambda",    "H3K4me3",    "negative binomial",        "updown",
#   "HMCan (H3K36me3)\nMACS (H3K4me3)",            "updown_poisson_constant_lambda",  "H3K36me3",   "HMCan (H3K36me3)\nMACS (H3K4me3)", "heuristics",
#   "HMCan (H3K36me3)\nMACS (H3K4me3)",            "updown_poisson_linear_lambda",    "H3K4me3",    "HMCan (H3K36me3)\nMACS (H3K4me3)", "heuristics"
# )

updown_poisson_vs_other <- tribble(
  ~target,                                        ~ref,                              ~type.peak,   ~lab2,                      ~lab,
  "PDPA_gaussian_constant_lambda",                "updown_poisson_constant_lambda",  "H3K36me3",   "gaussian",                 "unconstrained",
  "PDPA_poisson_constant_lambda",                 "updown_poisson_constant_lambda",  "H3K36me3",   "poisson",                  "unconstrained",
  "PDPA_negbin_constant_lambda_constant_phi",     "updown_poisson_constant_lambda",  "H3K36me3",   "negative binomial",        "unconstrained",
  "updown_gaussian_constant_lambda",              "updown_poisson_constant_lambda",  "H3K36me3",   "gaussian",                 "updown",
  "updown_negbin_constant_lambda_constant_phi",   "updown_poisson_constant_lambda",  "H3K36me3",   "negative binomial",        "updown",
  "PDPA_gaussian_linear_lambda",                  "updown_poisson_linear_lambda",    "H3K4me3",    "gaussian",                 "unconstrained",
  "PDPA_poisson_linear_lambda",                   "updown_poisson_linear_lambda",    "H3K4me3",    "poisson",                  "unconstrained",
  "PDPA_negbin_linear_lambda_constant_phi",       "updown_poisson_linear_lambda",    "H3K4me3",    "negative binomial",        "unconstrained",
  "updown_gaussian_linear_lambda",                "updown_poisson_linear_lambda",    "H3K4me3",    "gaussian",                 "updown",
  "updown_negbin_linear_lambda_constant_phi",     "updown_poisson_linear_lambda",    "H3K4me3",    "negative binomial",        "updown",
  "HMCan (H3K36me3)\nMACS (H3K4me3)",            "updown_poisson_constant_lambda",  "H3K36me3",   "HMCan (H3K36me3)\nMACS (H3K4me3)", "heuristics",
  "HMCan (H3K36me3)\nMACS (H3K4me3)",            "updown_poisson_linear_lambda",    "H3K4me3",    "HMCan (H3K36me3)\nMACS (H3K4me3)", "heuristics"
)

res_cv_all_updown_poisson_vs_other <- do.call(rbind,pmap(updown_poisson_vs_other, function(target, ref, type.peak, lab, lab2){
  type.peak_ <- type.peak
  cbind(
    res_cv_all[model == target & type.peak == type.peak_], 
    diff.accuracy = res_cv_all[model == target & type.peak == type.peak_]$accuracy - res_cv_all[model == ref & type.peak == type.peak_]$accuracy,
    lab = lab,
    lab2 = lab2
  )
}))

# target =   "PDPA_gaussian_linear_lambda"
# ref =  "updown_poisson_linear_lambda"
# type.peak = "H3K4me3"
# lab = "gaussian"
# lab2 = "PDPA"
target =   "PDPA_gaussian_linear_lambda"
ref =  "updown_poisson_linear_lambda"
type.peak = "H3K4me3"
lab = "gaussian"
lab2 = "unconstrained"

signi_updown_poisson_vs_other <- do.call(rbind,pmap(updown_poisson_vs_other, function(target, ref, type.peak, lab, lab2){
  type.peak_ <- type.peak
  r <- res_cv_all[(model == target | model == ref) & type.peak==type.peak_,]
  label <- mean_cl_normal(res_cv_all[model == target & type.peak==type.peak_]$accuracy - res_cv_all[model == ref & type.peak==type.peak_]$accuracy)
  res <- t.test(res_cv_all[model == target & type.peak==type.peak_]$accuracy, res_cv_all[model == ref & type.peak==type.peak_]$accuracy, paired=T)
  p.value <- res$p.value
  data.table(
    diff.accuracy = NaN,
    type.peak = type.peak,
    lab = lab,
    lab2 = lab2,
    p.value = p.value,
    x = label$y,
    xmin = label$ymin,
    xmax = label$ymax
  )
}))

signi_updown_poisson_vs_other[, adj.p.value := p.adjust(p.value, method="fdr")]
signi_updown_poisson_vs_other[,label := signi(adj.p.value)]

signi_updown_poisson_vs_other[,type.peak := ifelse(type.peak=="H3K36me3","H3K36me3 (broad peaks)","H3K4me3 (sharp peaks)" )]
res_cv_all_updown_poisson_vs_other[,type.peak := ifelse(type.peak=="H3K36me3","H3K36me3 (broad peaks)","H3K4me3 (sharp peaks)" )]

g <- ggplot(aes(y=as.factor(lab2), x=diff.accuracy), data=res_cv_all_updown_poisson_vs_other)+
  theme_bw()+
  facet_grid(lab~type.peak, scales ="free")+
  geom_vline(xintercept = 0, color="Royalblue")+
  geom_errorbar(aes(xmin = xmin, xmax=xmax), data= signi_updown_poisson_vs_other, color="red",width=0.5, size=1)+
  geom_point(aes(x = x), data= signi_updown_poisson_vs_other, color="red")+
  geom_text(aes(x = xmin -1, y=as.factor(lab2), label=label), data= signi_updown_poisson_vs_other)+
  xlab("accuracy(target model)% - accuracy(reference model)% on test folds")+
  ylab("target model")+
  theme(
    strip.background = element_rect(fill="grey95"), 
    strip.text.y = element_text(angle = 0),
    panel.margin=grid::unit(0, "lines"),
    text = element_text(size=15)
  )+
  ggtitle("reference model : updown poisson")

pdf(file = "figures/figure_updown_vs_other.pdf",
  width = 8,
  height = 3,
)
print(g)
dev.off()
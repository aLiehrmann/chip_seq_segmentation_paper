library(data.table)
library(purrr)
library(furrr)
library(ggplot2)
library(stringr)
library(aricode)

load(file="data/res_ari_heuristics.RData")
load(file="data/res_nid_heuristics.RData")
load(file = "data/res_ari.RData")
load(file = "data/res_nid.RData")


setkey(res, chunk.name, cell.type)
setkey(res_heuristics, chunk.name, cell.type)
res_all <- res[res_heuristics]

res_ <- rbind(res, res_heuristics)

g2 <- ggplot(data=res_, aes(x=ari, factor(model))) + 
facet_grid(.~type_peak, scales="free")+
geom_boxplot()+
theme_bw()+
xlab("NID")+
  ylab("model")+
  theme(
    strip.background = element_rect(fill="grey95"), 
    strip.text.y = element_text(angle = 0),
    panel.margin=grid::unit(0, "lines"),
    text = element_text(size=15)
  )

ggsave(g2, filename="figures/comparisons_robustness_models_boxplot.pdf", width=10, height=4)


res_all <- res_all[, list(diff.ari = (ari-i.ari), model=model, typePeak=type_peak, ari, i.ari)]

signi_diff_ari <- res_all[, mean_cl_normal(diff.ari),by=list(model, typePeak)]

names(signi_diff_ari) <- c("model","typePeak","x","xmin","xmax")
signi_diff_ari <- cbind(signi_diff_ari, p.value=res_all[,list(p.value=t.test(ari,  i.ari, paired=T)$p.value),by=list(model, typePeak)]$p.value)
signi_diff_ari[, adj.p.value := p.adjust(p.value, method="fdr")]
signi_diff_ari[,label := signi(adj.p.value)]

signi_diff_ari[,noise:=str_extract(model,"(?<= ).*")]
signi_diff_ari[,model:=str_extract(model,".*?(?= )")]

g <- ggplot( data= signi_diff_ari, aes(y=factor(noise),x=x, xmin=xmin, xmax=xmax))+
facet_grid(model~typePeak, scales="free")+
  geom_vline(xintercept = 0, color="Royalblue")+
  geom_errorbar(color="red",width=0.5, size=1)+
  geom_point(color="red")+
  geom_text(aes(x = xmin -0.03, label=label))+
  theme_bw()+
  theme(
    strip.background = element_rect(fill="grey95"), 
    strip.text.y = element_text(angle = 0),
    panel.margin=grid::unit(0, "lines"),
    text = element_text(size=15)
  )+
  xlab("NID(target model) - NID(reference model)")+
  ylab("target model")+
  ggtitle("reference model : HMCan (H3K36me3) & MACS (H3K4me3) heuristics")


pdf(file = "figures/comparison_robustness_models.pdf",
  width = 9.5,
  height = 3,
)
print(g)
dev.off()
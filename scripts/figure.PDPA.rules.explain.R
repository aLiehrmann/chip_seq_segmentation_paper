library(data.table)
library(ggplot2)

j = 6
s = 6
f.counts.H3K4me3 <- "data/chunks/H3K36me3_AM_immune/12/counts.RData"
load(f.counts.H3K4me3)
sample.id = s
df.counts.H3K4me3 <- split(counts, counts$sample.id)[[sample.id]]

f.regions <- "data/chunks/H3K36me3_AM_immune/12/regions.RData"
load(f.regions)
(df.regions <- split(regions, regions$sample.id)[[sample.id]])

f.models.H3K4me3.PDPA <- "data/chunks/H3K36me3_AM_immune/12/PDPA_poisson.RData"
sample.name <- df.counts.H3K4me3$sample.id[1]

load(f.models.H3K4me3.PDPA)
models.H3K4me3.PDPA <- dt.models[sample==sample.name]

#8
fit <- models.H3K4me3.PDPA[j,]$fit[[1]]
cp <- df.counts.H3K4me3$chromEnd[c(1,fit$changepoints)]
cp <- cp[-(1:14)]
n.cp <- length(cp)
m <- fit$parameters[-(1:14)]

df.cp <- data.table(x = cp[-c(1,n.cp)])
df.m <- data.table(
  x = cp[-n.cp],
  xend = cp[-1],
  y = m,
  yend = m
)

df.segments <- data.table(
  xmin = c(cp[c(3,2,2)]),
  xmax = c(cp[c(4,5,4)]),
  ymin= c(-3.5, -8.5, -13.5),
  ymax= c(-3.5, -8.5, -13.5)
)


g <- ggplot() +
  geom_path(
    aes(x=chromEnd/10e3, y=coverage),
    data = df.counts.H3K4me3,
    color="grey80"
  ) +
  geom_hline(
    yintercept=c(-0.5,-6,-11),
    color="grey70",
  )+
  geom_segment(
    aes(x=x/10e3,xend=xend/10e3,y=y,yend=yend),
    data = df.m,
    size=2,
    color="royalblue2"
  ) +
  geom_segment(
    aes(x=x/10e3, xend=x/10e3, y=-Inf, yend=55),
    data = df.cp,
    size = 1.5,
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    "text",
    x = df.cp$x[1]/10e3,
    y = 60,
    color="red",
    label = "Start1",
    size=6
  )+
  annotate(
    "text",
    x = df.cp$x[2]/10e3,
    y = 60,
    color="red",
    label = "Start2",
    size=6
  )+
  annotate(
    "text",
    x = df.cp$x[3]/10e3,
    y = 60,
    color="red",
    label = "End1",
    size=6
  )+
  annotate(
    "text",
    x = df.cp$x[4]/10e3,
    y = 60,
    color="red",
    label = "End2",
    size=6
  )+
  geom_segment(
    aes(
      x=xmin/10e3,
      y=ymin,
      xend=xmax/10e3,
      yend=ymax),
    lineend = "round",
    data = df.segments,
    size=2
  )+
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = 6282,
      ymin = -Inf,
      ymax = -0.5,
    ),
    fill = "grey95",
    color="black"
  )+
  geom_text(aes(
    x = 6282-0.6,
    y = -8.5,
    label = "rules"),
    angle=90,
    size=6
  )+
  geom_text(aes(
    x = 6282+2,
    y = -3.5,
    label = "(thinnest peak)"),
    size=4
  )+
  geom_text(aes(
    x = 6282+2,
    y = -8.5,
    label = "(largest peak)"),
    size=4
  )+
  geom_text(aes(
    x = 6282+2,
    y = -13.5,
    label = "(max jump)"),
    size=4
  )+
  theme_bw()+
  theme(
    text = element_text(
    size=18)
  )+
  ylab("aligned sequence\nreads")+
  xlab("position on chromosome (kb : kilo bases)")+
  coord_cartesian(xlim = c(
    6282,
    6305
  )) 

pdf("figures/figure_PDPA_rules_explain.pdf",
  width = 9.5,
  height = 4,
)
print(g)
dev.off()
library(data.table)
library(ggplot2)
library(ggnewscale)


f.counts <- "data/chunks/H3K36me3_AM_immune/13/counts.RData"
f.regions <- "data/chunks/H3K36me3_AM_immune/13/regions.RData"
f.models <- "data/chunks/H3K36me3_AM_immune/13/updown_poisson.RData"
load(f.counts)
load(f.regions)
load(f.models)
sample.id = 10
df.counts <- split(counts, counts$sample.id)[[sample.id]]
sample.name <- df.counts$sample.id[1]
df.regions <- split(regions, regions$sample.id)[[sample.id]]
models <- dt.models[sample==sample.name]
df.regions

# CP (POINTS)

df.cp.m1 <- data.table(
  cp.x= df.counts$chromEnd[
    models[2,]$fit[[1]]$changepoints[-3]
  ]/10e3,
  cp.y = -3.5
)

df.cp.m2 <- data.table(
  cp.x = df.counts$chromEnd[
    models[3,]$fit[[1]]$changepoints[-5]
  ]/10e3,
  cp.y = -8.5
)

df.cp.m3 <- data.table(
  cp.x = df.counts$chromEnd[
    models[4,]$fit[[1]]$changepoints[-7]
  ]/10e3,
  cp.y = -13.5
)
df.cp.m3[3,1] <- df.cp.m3[3,1]-0.1
df.cp.m3[4,1] <- df.cp.m3[4,1]+0.1
df.cp <- rbind(
  df.cp.m1,
  df.cp.m2,
  df.cp.m3
)


# CP (SEGMENTS)


df.segments.m1 <- data.table(
  xmin = df.counts$chromEnd[
    models[2,]$fit[[1]]$changepoints[1]
  ]/10e3,
  xmax = df.counts$chromEnd[
    models[2,]$fit[[1]]$changepoints[2]
  ]/10e3,
  ymin=-3.5,
  ymax=-3.5
)

df.segments.m2 <- data.table(
  xmin = df.counts$chromEnd[
    models[3,]$fit[[1]]$changepoints[c(1,3)]
  ]/10e3,
  xmax = df.counts$chromEnd[
    models[3,]$fit[[1]]$changepoints[c(2,4)]
  ]/10e3,
  ymin=-8.5,
  ymax=-8.5
)

df.segments.m3 <- data.table(
  xmin = df.counts$chromEnd[
    models[4,]$fit[[1]]$changepoints[c(1,3,5)]
  ]/10e3,
  xmax = df.counts$chromEnd[
    models[4,]$fit[[1]]$changepoints[c(2,4,6)]
  ]/10e3,
  ymin=-13.5,
  ymax=-13.5
)

df.segments <- rbind(
  df.segments.m1,
  df.segments.m2,
  df.segments.m3
)

# STATUS 

df.status.m1 <- data.table(
  xmin = df.regions$chromStart/10e3,
  xmax = df.regions$chromEnd/10e3, 
  ymin = -6,
  ymax = -1,
  status = c(
    "correct",
    "correct",
    "FN",
    "FN",
    "correct")
)

df.status.m2 <- data.table(
  xmin = df.regions$chromStart/10e3,
  xmax = df.regions$chromEnd/10e3, 
  ymin = -11,
  ymax = -6,
  status = c(
    "correct",
    "correct",
    "correct",
    "correct",
    "correct")
)


df.status.m3 <- data.table(
  xmin = df.regions$chromStart/10e3,
  xmax = df.regions$chromEnd/10e3, 
  ymin = -16,
  ymax = -11,
  status = c(
    "correct",
    "correct",
    "FP",
    "correct",
    "correct")
)

df.status <- rbind(
  df.status.m1,
  df.status.m2,
  df.status.m3
)
# PLOT

ann.colors <- c(
  noPeaks="grey90",
  peakEnd="royalblue",
  peakStart="steelblue1",
  peaks="#a445ee"
)

status.colors <- c(
  correct="forestgreen",
  FP="red2",
  FN="orange")

g <- ggplot()+
  geom_path(
    aes(x=chromEnd/10e3, y=coverage), 
    data=df.counts,
    size=0.6,
    color="grey40"
  )+
  geom_rect(
    aes(
      xmin=chromStart/10e3,
      xmax=chromEnd/10e3, 
      ymin=0, 
      ymax=Inf, 
      fill = annotation
    ),
    data = df.regions,
    alpha=0.7
  )+
  scale_fill_manual(values=ann.colors)+
  new_scale_fill()+
  geom_vline(
    aes(xintercept=c(
      df.regions$chromStart/10e3, 
      df.regions$chromEnd/10e3
      )
    ),
    color="grey70"
  )+
  scale_y_continuous(breaks=c(0,10,20,30,40))+
  geom_hline(
    yintercept=c(-1,-6,-11),
    color="grey70"
  )+
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = min(df.counts$chromEnd/10e3),
      ymin = -Inf,
      ymax = -1,
    ),
    fill = "grey95",
    color="black"
  )+
  geom_rect(
    aes(
      xmin=xmin,
      ymin=ymin,
      xmax=xmax,
      ymax=ymax,
      fill = status
    ),
    data=df.status,
    alpha=0.7
  )+
  scale_fill_manual(values=status.colors)+
  geom_segment(
    aes(
      x=xmin,
      y=ymin,
      xend=xmax,
      yend=ymax),
    lineend = "round",
    data = df.segments,
    size=2
  )+
  geom_text(aes(
    x = min(df.counts$chromEnd/10e3)-0.5,
    y = -8.5,
    label = "models"),
    angle=90,
    size=6
  )+
  geom_text(aes(
    x = min(df.counts$chromEnd/10e3)+1,
    y = -3.5,
    label = "(1 peak)"),
    size=5
  )+
  geom_text(aes(
    x = min(df.counts$chromEnd/10e3)+1,
    y = -8.5,
    label = "(2 peaks)"),
    size=5
  )+
  geom_text(aes(
    x = min(df.counts$chromEnd/10e3)+1,
    y = -13.5,
    label = "(3 peaks)"),
    size=5
  )+
  theme_bw()+
  theme(
    text = element_text(
    size=18)
  )+
  ylab("aligned sequence\nreads")+
  xlab("position on chromosome (kb : kilo bases)")+
  coord_cartesian(
    ylim= c(-14.5, max(df.counts$coverage))+1,
    xlim = c(
      min(df.counts$chromEnd/10e3)-0.22, 
      max(df.counts$chromEnd/10e3)-0.3
    )
  )

pdf("figures/figure_labeled_data.pdf",
  width = 11 ,
  height = 4
)
print(g)
dev.off()


# library(tikzDevice)
# tikz("figures/FIG_labeled_data.tex",
#   width = 11 ,
#   height = 4
# )
# print(g)
# dev.off()
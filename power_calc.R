library(ggplot2)
library(gridExtra)
library(grid)


source("functions.R")
source("avengeme.R")
a <- expand.grid(
    vg1=c(0.01, 0.05, 0.1, 0.15),  
    prop=seq(0, 1, 0.01),
    prev=c(0.07, 0.12)
)
a$cov12 <- a$vg1 * a$prop



for(i in 1:nrow(a))
{
    a$power[i] <- polygenescore(
        nsnp = 100000,  
        n = c(100000,8000),  
        vg1 = a$vg1[i],
        cov12 = a$cov12[i],
        binary = c(FALSE, TRUE),
        prevalence = c(0, a$prev[i]),
        sampling = c(0, a$prev[i])
        )$power
}


a$lab <- as.character(a$prev)
a$lab[a$lab == "0.07"] <- "SA"
a$lab[a$lab == "0.12"] <- "NSSI"

p1 <- ggplot(a, aes(x=prop, y=power)) +
geom_line(aes(colour=as.factor(vg1))) +
facet_grid(. ~ lab) +
labs(x="Hypothesised genetic correlation", colour="Variance explained\nin trait 1")
ggsave("pgs_power.pdf", height=4, width=6)


b <- expand.grid(
    h21=c(0.1, 0.3, 0.5, 0.7),  
    prop=seq(0, 1, 0.01),
    prev=c(0.07, 0.12)
)

for(i in 1:nrow(b))
{
    b$power[i] <- calcBiQtCc(
        n=5000,
        ncase=5000 * b$prev[i],
        ncontrol=5000 * (1 - b$prop[i]),
        hsq1=b$h21[i],
        hsq2=0.5,
        K=b$prop[i],
        rg=b$prop[i],
        overlap=TRUE
    )$power
}

b$lab <- as.character(b$prev)
b$lab[b$lab == "0.07"] <- "SA"
b$lab[b$lab == "0.12"] <- "NSSI"


p2 <- ggplot(b, aes(x=prop, y=power)) +
geom_line(aes(colour=as.factor(h21))) +
facet_grid(. ~ lab) +
labs(x="Hypothesised genetic correlation", colour="Heritability of trait 1")
ggsave("gcta_power.pdf", height=4, width=6)


pdf("graphs.pdf", width=6, height=4)
grid.arrange(p1 + labs(y="PGS power", x=NULL), p2 + labs(y="GREML power"), ncol=1)
dev.off()

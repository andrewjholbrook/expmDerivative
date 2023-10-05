setwd("~/expmDerivative/")
library(ggplot2)

exact <- readRDS("data/lowDimSimExact.rds")
approx <- readRDS("data/lowDimSimApprox.rds")
correct <- readRDS("data/lowDimSimCorrected.rds")

plot(exact[[2]],type="l")
lines(approx[[2]],col="red")
lines(correct[[2]],col="green")

# flatten lists
df_ex <- matrix(0,100000,25)
for (i in 1:100000) {
  exact[[1]][[i]] <- exp(exact[[1]][[i]])
  for(j in 1:5) {
    exact[[1]][[i]][j,j] <- -(exact[[1]][[i]][j,j] - sum(exact[[1]][[i]][j,]))
  }
  df_ex[i,] <- as.vector(exact[[1]][[i]])
}
df_ex <- data.frame(df_ex)
df_ex$method <- "exact"

df_ap <- matrix(0,100000,25)
for (i in 1:100000) {
  approx[[1]][[i]] <- exp(approx[[1]][[i]])
  for(j in 1:5) {
    approx[[1]][[i]][j,j] <- -(approx[[1]][[i]][j,j] - sum(approx[[1]][[i]][j,]))
  }
  df_ap[i,] <- as.vector(approx[[1]][[i]])
}
df_ap <- data.frame(df_ap)
df_ap$method <- "approx"

df_co <- matrix(0,100000,25)
for (i in 1:100000) {
  correct[[1]][[i]] <- exp(correct[[1]][[i]])
  for(j in 1:5) {
    correct[[1]][[i]][j,j] <- -(correct[[1]][[i]][j,j] - sum(correct[[1]][[i]][j,]))
  }
  df_co[i,] <- as.vector(correct[[1]][[i]])
}
df_co <- data.frame(df_co)
df_co$method <- "correct"

df <- rbind(df_ap,df_ex,df_co)

df$Gradient <- factor(df$method)
levels(df$Gradient) <- c("Naive","Corrected","Exact")

df <- reshape2::melt(df,value.name = "value",measure.vars=1:25)
df$variable <- factor(df$variable)
levels(df$variable) <- c("1:1","1:2","1:3","1:4","1:5",
                         "2:1","2:2","2:3","2:4","2:5",
                         "3:1","3:2","3:3","3:4","3:5",
                         "4:1","4:2","4:3","4:4","4:5",
                         "5:1","5:2","5:3","5:4","5:5")

set.seed(1)
S <- 5
maxIts <- 100000
initialStates <- sample(x=1:S,size=20,replace=TRUE) # only doing this to keep seed the same
Q       <- matrix(exp(rnorm(S^2)),S,S)
diag(Q) <- 0
diag(Q) <- - rowSums(Q)
Q <- as.vector(Q)
df$Truth <- abs(rep(Q,each=300000))

gg <- ggplot(df,aes(x=value,fill=Gradient,color=Gradient)) +
  # geom_errorbarh(aes(color=Gradient,y=0.4,xmin=stat_quantile(quantiles = 0.025),
  #                    xmax=stat_quantile(quantiles = 0.975))) +
  geom_density(alpha=0.3,adjust=1) +
  geom_vline(aes(xintercept = Truth),color="red") +
  facet_wrap(~variable) +
  scale_fill_manual(values = ggthemes::stata_pal("s2color")(15)[c(3,4,7)]) +
  scale_color_manual(values = ggthemes::stata_pal("s2color")(15)[c(3,4,7)]) +
  #ylim(c(0,0.5)) +
  coord_cartesian(ylim=c(0,0.5),xlim=c(0,25)) +
  xlab("Value") + ylab("Density") + #+ ggtitle("Empirical posterior distributions") +
  theme_minimal()
gg

ggsave(filename="figures/lowDimVis.pdf",plot=gg,device="pdf",width = 7,height=6)

system2(command = "pdfcrop",
        args    = c("~/expmDerivative/figures/lowDimVis.pdf",
                    "~/expmDerivative/figures/lowDimVis.pdf")
)
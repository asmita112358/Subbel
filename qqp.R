##code for plotting qqplots
library(ggplot2)
library(data.table)
library(ggpubr)
library(Haplin)
c1 = read.csv("size_c1.csv")[,-1]
c2 = read.csv("size_c2.csv")[,-1]
c3 = read.csv("size_c3.csv")[,-1]
c4 = read.csv("size_c4.csv")[,-1]
N_rep = nrow(c1)
omega = 0.05
allpvals = list(c1,c2,c3,c4)

qqp = list()
barp = list()
title = rep(c("Sparse null", "Dense null"), each = 2)
size <- se_size <- tmp.data <- list()

#yup <- c(0.2, 0.25, 0.2, 0.25)

for(case in seq_along(allpvals)){
  x_loc <- -log10(0.05)
  y_loc <- -log10(0.05)
  mc = allpvals[[case]]
  N_rep = nrow(mc)
  log.pvals = -log10(apply(mc, 2, sort))
  cat(max(log.pvals))
  log.theo = -log10((1:N_rep)/(N_rep+1))
  data = data.frame(log.theo,log.pvals)
  colnames(data) <- c("theo_q", "Subbel", "ABtest", "MaxP", "Sobel")
  data.plot <- data.table(data)
  .lim <- c(0, max(log.pvals)) * 1.05
  qqp[[case]] <- data.plot %>%
    melt(id.vars = 1,
         variable.name = "Method") %>%
    ggplot(aes(x = theo_q, y = value,
               group = Method, color = Method,
               shape = Method)) +
    geom_point() +
    geom_abline() +
    geom_segment(aes(x = x_loc, y = 0, xend = x_loc, yend = x_loc),
                 linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, y = y_loc, xend = y_loc, yend = y_loc),
                 linetype = "dotted", color = "black", linewidth = 0.5) +
    labs(title = title[case],
         x = "Expected P-value (-log10 scale)", y = "") +
    scale_x_continuous(limits = c(0,4), breaks = 0:4, expand = c(0,0)) +
    scale_y_continuous(limits = c(0,4), breaks = 0:4, expand = c(0,0)) +
    theme_minimal()+
    theme(panel.border = element_rect(color = "black", fill = NA))
  
  qqp[[case]]
}

qqp_low <- ggarrange(qqp[[1]], qqp[[3]], nrow = 1,align = "h", common.legend = T)
final_qqp <- annotate_figure(qqp_low, left = text_grob("Observed P-value (-log10 scale)", rot = 90, vjust = 1, size = 14))
final_qqp

ggsave("qqp_r1.png", final_qqp, device = "png",
       width = 8, height = 4.5, units = "in", bg = "white",
       dpi = 1200)

qqp_high <- ggarrange(qqp[[2]], qqp[[4]], nrow = 1,align = "h", common.legend = T)
final_qqp <- annotate_figure(qqp_high, left = text_grob("Observed P-value (-log10 scale)", rot = 90, vjust = 1, size = 14))
final_qqp

ggsave("qqp_r2.png", final_qqp, device = "png",
       width = 8, height = 4.5, units = "in", bg = "white",
       dpi = 1200)

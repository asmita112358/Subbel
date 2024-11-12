library(data.table)
library(ggplot2)
library(ggpubr)
library(Haplin)
library(ggpattern)
##Subbel size plots

rm(list = ls())
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

yup <- c(0.2, 0.25, 0.2, 0.25)

for(case in seq_along(allpvals)){
  
  mc = allpvals[[case]]
  
  theo_q = seq(0.01,1, 0.01)
  samp_q = matrix(nrow = 100, ncol = 4)
  for(i in 1:4)
  {
    samp_q[,i] = quantile(mc[,i], probs = theo_q)
  }
  data = data.frame(-log10(theo_q), -log10(samp_q))
  colnames(data) <- c("theo_q", "Subbel", "ABtest", "MaxP", "Sobel")
  data.plot <- data.table(data)
  
  size[[case]] = apply((mc<=omega), 2, mean)
  se_size[[case]] = sqrt(size[[case]]*(1 - size[[case]])/ N_rep)
  
  tmp.data[[case]] <- data.table(Method = factor(colnames(data.plot)[2:5],
                                                 levels = c("Subbel", "ABtest",
                                                            "MaxP", "Sobel")),
                                 Size = size[[case]],
                                 SE_size = se_size[[case]])
  
  barp[[case]] <-  tmp.data[[case]] %>%
    ggplot(aes(x = Method, y = Size, fill = Method)) +
    geom_bar(stat = "identity") +
    scale_pattern_manual(
      "Method",
      values = c("none", "stripe", "none", "none")
    )+
    geom_col_pattern(
      aes(pattern = Method),
      position = "dodge",
      pattern_angle = 45,
      pattern_density = .1,
      pattern_spacing = .04,
      pattern_fill = "black",
      color = "black"
    ) +
    geom_errorbar(aes(ymin = Size - 1.96*SE_size, 
                      ymax = Size + 1.96*SE_size),
                  width = 0.2,position = position_dodge(width = 0.9)) + 
    geom_hline(yintercept = 0.05) +
    labs(title = title[case]) +
    ylim(-0.01, yup[case]) +
    xlab("")+
    ylab("")+
    theme_minimal()+
    theme(panel.border = element_rect(color = "black", fill = NA))
  
}

barp_low <- ggarrange(barp[[1]], barp[[3]], nrow = 1,align = "h", common.legend = T)
final_barp <- annotate_figure(barp_low, left = text_grob("Size", rot = 90, vjust = 1, size = 14))
final_barp

ggsave("barp_r1.png", final_barp, device = "png",
       width = 8, height = 4.5, units = "in", bg = "white",
       dpi = 1200)
barp_high <- ggarrange(barp[[2]], barp[[4]], nrow = 1,align = "h", common.legend = T)
final_barp <- annotate_figure(barp_high, left = text_grob("Size", rot = 90, vjust = 1, size = 14))
final_barp

ggsave("barp_r2.png", final_barp, device = "png",
       width = 8, height = 4.5, units = "in", bg = "white",
       dpi = 1200)

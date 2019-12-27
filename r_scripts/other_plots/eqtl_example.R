suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
options(bitmapType = "cairo")
set.seed(123)


data1 <- data.frame(x = c(rep("XX", 10), rep("XY", 10), rep("YY", 10)), y = c(c(1:10), c(6:15), c(11:20))) %>%
  as.tibble()

plot1 <- ggplot(data = data1, mapping=aes(x=x, y=y, fill=x)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size=35),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    aspect.ratio=1) +
  scale_fill_brewer(type = "seq", palette = "Greys") +
labs(x = "Genotypes", y = "Expression")
ggsave(filename="eqtl_example.png", plot=plot1, height=3, width=3, path = "./plots/other_plots/")
ggsave(filename="eqtl_example.pdf", plot=plot1, height=3, width=3, path = "./plots/other_plots/")
unlink("eqtl_example.png")
unlink("eqtl_example.pdf")


data2 <- data.frame(x1 = c(1:2500), x2 = c(rep("1", 500), rep("2", 500), rep("3", 500), rep("4", 500), rep("5", 500)), y = runif(2500, min = 0.05, max = 1)) %>%
  as.tibble() %>%
  mutate(y = replace(y, sample(1:2500, 10), runif(10, min = 0.01, max = 0.05))) %>%
  mutate(y = -log10(y))


plot2 <- ggplot(data = data2, mapping=aes(x=x1, y=y, color=x2)) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), colour="black", linetype=2, size = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size=22),
    axis.title.y = element_text(size=22, margin = margin(r = 10)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=16, color="black"),
    plot.title = element_text(size=22),
    legend.position = "none",
    aspect.ratio = 1) +
    scale_y_continuous(breaks = c(-log10(0.05)), labels = c("0.05")) +
labs(x = "Chromosome", y = "P-value", title="GWAS")
ggsave(filename="gwas_example1.png", plot=plot2, height=3, width=3, path = "./plots/other_plots/")
unlink("gwas_example1.png")



data3 <- data.frame(x1 = c(1:2500), x2 = c(rep("1", 500), rep("2", 500), rep("3", 500), rep("4", 500), rep("5", 500)), y = runif(2500, min = 0.05, max = 1)) %>%
  as.tibble() %>%
  mutate(y = replace(y, sample(1:2500, 2), runif(2, min = 0.01, max = 0.05))) %>%
  mutate(y = -log10(y))


plot3 <- ggplot(data = data3, mapping=aes(x=x1, y=y, color=x2)) +
  geom_point() +
  geom_hline(yintercept=-log10(0.05), colour="black", linetype=2, size = 0.5) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size=22),
    axis.title.y = element_text(size=22, margin = margin(r = 10)),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=16, color="black"),
    plot.title = element_text(size=22),
    legend.position = "none",
    aspect.ratio = 1) +
    scale_y_continuous(breaks = c(-log10(0.05)), labels = c("0.05")) +
labs(x = "Chromosome", y = "P-value", title="GWAS")
ggsave(filename="gwas_example2.png", plot=plot3, height=3, width=3, path = "./plots/other_plots/")
unlink("gwas_example2.png")

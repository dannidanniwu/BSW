dp <- dd[, .(avg = mean(y)), keyby = .(A, site, k)]

ggplot(data = dp, aes(x = k, y = avg)) +
  geom_line(aes(group = site, color = factor(A))) +
  scale_color_manual(values = c("#a53e3e", "#3ea5a5"),
                     guide = guide_legend(reverse = TRUE),
                     labels = c("no treatment", "treatment")) +
  ylab("average Y") +
  xlab("period (k)") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank())

ggplot(data = dd, aes(x = k, y = y)) +
  geom_point(aes(color = factor(A, labels = c("Control", "Intervention"))),
             size = 0.1) +
  scale_color_manual(values = c("#d07b7c", "#7ba7d0")) +
  facet_wrap(~site, ncol = 8) +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  xlim(c(0, 24)) +
  guides(color = guide_legend(override.aes = list(size = 2)))
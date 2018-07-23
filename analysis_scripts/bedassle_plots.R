library(BEDASSLE)
library(cowplot)

# plots posterior distribtions
# lot of steps in here are annoying and manual, these aren't ideal plotting methods!

# all pops
load("bedassle_output/bed_tp_run10_main_MCMC_output1.Robj")
# aD should be multiplied by 164.7215, the standard deviation of geo dist
# aE2 should be multiplied by 7.817393 the standard deviation of ppt diff

ae1_ad = aE[1,]/(aD*164.7215*100)

# ae1_ad_trim = ae1_ad[100:10000]
ae1_ad_trim = ae1_ad[2001:10000]

ae2_ad = (aE[2,]*7.817393*10)/(aD*164.7215*100)

# ae2_ad_trim = ae2_ad[100:10000]
ae2_ad_trim = ae2_ad[2001:10000]

load("bedassle_output/bed_tp_run12_north_MCMC_output1.Robj")
# aD should be multiplied by 164.7215
# aE2 should be multiplied by 7.817393

ae1_ad_n = aE[1,]/(aD*164.7215*100)

# ae1_ad_trim_n = ae1_ad_n[100:10000]
ae1_ad_trim_n = ae1_ad_n[2001:10000]

ae2_ad_n = (aE[2,]*7.817393*10)/(aD*164.7215*100)

# ae2_ad_trim_n = ae2_ad_n[100:10000]
ae2_ad_trim_n = ae2_ad_n[2001:10000]

load("bedassle_output/bed_tp_run14_center_MCMC_output1.Robj")
# aD should be multiplied by  164.7215
# aE2 should be multiplied by 7.817393

ae1_ad_c = aE[1,]/(aD*164.7215*100)

# ae1_ad_trim_c = ae1_ad_c[100:10000]
ae1_ad_trim_c = ae1_ad_c[2001:10000]


ae2_ad_c = (aE[2,]*7.817393*10)/(aD*164.7215*100)

# ae2_ad_trim_c = ae2_ad_c[100:10000]
ae2_ad_trim_c = ae2_ad_c[2001:10000]

main_ests = data.frame(ae1_ad_trim, ae2_ad_trim, ae1_ad_trim_n, ae2_ad_trim_n, ae1_ad_trim_c, ae2_ad_trim_c)

up1 = stats::quantile(ae1_ad_trim, prob = c(0.025, 0.975), na.rm = TRUE)[2]
dwn1 = stats::quantile(ae1_ad_trim, prob = c(0.025, 0.975), na.rm = TRUE)[1]
med1 = median(ae1_ad_trim, na.rm = TRUE)
med1print = round(med1, 9)

med1print

m1 =
  ggplot(main_ests, aes(x = ae1_ad_trim, y = ..scaled..)) +
  geom_density(fill = "lightblue") +
  geom_segment(aes(y = 0, yend =  0.15, x = up1, xend = up1), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.16, x = dwn1, xend = dwn1), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  1, x = med1, xend = med1), size = 0.75) +  
  ylab("") +
  scale_x_continuous(breaks = c(1e-7, 1.5e-7, 2e-7), labels = c(expression("1.0 x 10"*{}^{-7}*""),
                                                                expression("1.5 x 10"*{}^{-7}*""),
                                                                expression("2.0 x 10"*{}^{-7}*""))) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab(expression(frac("Effect of 1°C temperature difference", "Effect of 100 km geographic distance"))) #+
  # annotate("text", x = med1 + med1*0.3, y = 0.945, label = paste("'median ='~1.18~'x'~10^-7"), parse = TRUE, size = 4)

up2 = stats::quantile(ae2_ad_trim, prob = c(0.025, 0.975), na.rm = TRUE)[2]
dwn2 = stats::quantile(ae2_ad_trim, prob = c(0.025, 0.975), na.rm = TRUE)[1]
med2 = median(ae2_ad_trim, na.rm = TRUE)
med2print = round(med2, 7)

m2 =
  ggplot(main_ests, aes(x = ae2_ad_trim, y = ..scaled..)) +
  geom_density(fill = "lightblue") +
  geom_segment(aes(y = 0, yend =  0.033, x = up2, xend = up2), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.65, x = dwn2, xend = dwn2), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.682, x = med2, xend = med2), size = 0.75) +  
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  xlab(expression(frac("Effect of 10 mm precipitation difference", "Effect of 100 km geographic distance"))) +
  # annotate("text", x = med2*3, y = 0.72, label = paste("'median ='~5.84~'x'~10^-7"), parse = TRUE, size = 4) +
  scale_x_continuous(breaks = c(0, 2e-6, 4e-6, 6e-6), labels = c("0", expression("2.0 x 10"*{}^{-6}*""), 
                                                                expression("4.0 x 10"*{}^{-6}*""),
                                                                expression("6.0 x 10"*{}^{-6}*""))) 
plot_grid(m1, m2, labels = c("A", "B")) 
ggsave("figs_for_paper/bedassle_main.pdf", height = 4.5, width = 8.5)


up3 = stats::quantile(ae1_ad_trim_n, prob = c(0.025, 0.975), na.rm = TRUE)[2]
dwn3 = stats::quantile(ae1_ad_trim_n, prob = c(0.025, 0.975), na.rm = TRUE)[1]
med3 = median(ae1_ad_trim_n, na.rm = TRUE)
med3print = round(med3, 7)

n1 =
  ggplot(main_ests, aes(x = ae1_ad_trim_n, y = ..scaled..)) +
  geom_density(fill = "lightblue") +
  geom_segment(aes(y = 0, yend =  0.19, x = up3, xend = up3), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.26, x = dwn3, xend = dwn3), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.99, x = med3, xend = med3), size = 0.75) +  
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  +
  # annotate("text", x = med3*1.8, y = 0.97, label = paste("'median ='~5.89~'x'~10^-8"), parse = TRUE, size = 5)  +
  xlab("")  +
    scale_x_continuous(breaks = c(0, 5e-8, 1e-7), labels = c(expression("0"), expression("5.0 x 10"*{}^{-8}*""), 
                                                                   expression("1.0 x 10"*{}^{-7}*""))) #+
  # xlab(expression(frac("Effect of 1°C temperature difference", "Effect of 100 km geographic distance")))

up4 = stats::quantile(ae2_ad_trim_n, prob = c(0.025, 0.975), na.rm = TRUE)[2]
dwn4 = stats::quantile(ae2_ad_trim_n, prob = c(0.025, 0.975), na.rm = TRUE)[1]
med4 = median(ae2_ad_trim_n, na.rm = TRUE)
med4print = round(med4, 7)

n2 =
  ggplot(main_ests, aes(x = ae2_ad_trim_n, y = ..scaled..)) +
  geom_density(fill = "lightblue") +
  geom_segment(aes(y = 0, yend =  0.17, x = up4, xend = up4), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.49, x = dwn4, xend = dwn4), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  1, x = med4, xend = med4), size = 0.75) +  
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())  +
  # annotate("text", x = med4*2.3, y = 0.9, label = paste("'median ='~9.73~'x'~10^-6"), parse = TRUE, size = 5) +
  xlab("")  +
  scale_x_continuous(breaks = c(0, 1e-5, 2e-5), labels = c(expression("0"), expression("1.0 x 10"*{}^{-5}*""), 
                                                           expression("2.0 x 10"*{}^{-5}*""))) #+
  # xlab(expression(frac("Effect of 10 mm precipitation difference", "Effect of 100 km geographic distance")))

up5 = stats::quantile(ae1_ad_trim_c, prob = c(0.025, 0.975), na.rm = TRUE)[2]
dwn5 = stats::quantile(ae1_ad_trim_c, prob = c(0.025, 0.975), na.rm = TRUE)[1]
med5 = median(ae1_ad_trim_c, na.rm = TRUE)
med5print = round(med5, 7)

c1 =
  ggplot(main_ests, aes(x = ae1_ad_trim_c, y = ..scaled..)) +
  geom_density(fill = "lightblue") +
  geom_segment(aes(y = 0, yend =  0.22, x = up5, xend = up5), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.38, x = dwn5, xend = dwn5), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  1, x = med5, xend = med5), size = 0.75) +  
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())   +
    scale_x_continuous(breaks = c(0, 2e-7, 4e-7), labels = c(expression("0"), expression("2.0 x 10"*{}^{-7}*""), 
                                                             expression("4.0 x 10"*{}^{-7}*""))) +
  # annotate("text", x = med5*1.8, y = 0.975, label = paste("'median ='~2.34~'x'~10^-7"), parse = TRUE, size = 5) +
  xlab(expression(frac("Effect of 1°C temperature difference", "Effect of 100 km geographic distance")))

up6 = stats::quantile(ae2_ad_trim_c, prob = c(0.025, 0.975), na.rm = TRUE)[2]
dwn6 = stats::quantile(ae2_ad_trim_c, prob = c(0.025, 0.975), na.rm = TRUE)[1]
med6 = median(ae2_ad_trim_c, na.rm = TRUE)
med6print = round(med6, 7)

c2 =
  ggplot(main_ests, aes(x = ae2_ad_trim_c, y = ..scaled..)) +
  geom_density(fill = "lightblue") +
  geom_segment(aes(y = 0, yend =  0.03, x = up6, xend = up6), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.65, x = dwn6, xend = dwn6), linetype = "dashed", size = 0.75) +  
  geom_segment(aes(y = 0, yend =  0.64, x = med6, xend = med6), size = 0.75) +  
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(breaks = c(0, 2e-6, 4e-6, 6e-6), labels = c(expression("0"), expression("2.0 x 10"*{}^{-6}*""), 
                                                           expression("4.0 x 10"*{}^{-6}*""), 
                                                           expression("6.0 x 10"*{}^{-6}*"")))  +
  # annotate("text", x = med6*2.9, y = 0.72, label = paste("'median ='~9.46~'x'~10^-7"), parse = TRUE, size = 5) +
  xlab(expression(frac("Effect of 10 mm precipitation difference", "Effect of 100 km geographic distance")))

sub1 = plot_grid(n1, n2, ncol = 2, labels = c("A", "B"))

sub2 = plot_grid(c1, c2, ncol = 2, labels = c("C", "D"))

n_title = ggdraw() + draw_label("Northern populations", fontface='bold', hjust = 1.9)
c_title = ggdraw() + draw_label("Central populations", fontface='bold', hjust = 2.05)

plot_grid(sub1, sub2, ncol = 1, rel_heights = c(0.9, 1))

plot_grid(n_title, sub1, c_title, sub2, ncol = 1, rel_heights = c(0.1, 0.9, 0.1, 1))
ggsave("figs_for_paper/bedassle_regional.pdf", width = 8.5, height = 8.5)  


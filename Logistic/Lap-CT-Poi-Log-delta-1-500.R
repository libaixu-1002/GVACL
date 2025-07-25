
library(ggplot2)
log10_labels <- function(x) {
  parse(text = paste0("10^", log10(x)))  # 转换为 10^x 格式
}
M <- c(18,20,30,50,80,100,150,200,300,400)
N 
#N <- c(12,20,30,50,80,100,150,200,300,400)
plot1 <- 
  ggplot()+                   # 添加连线
  geom_point(data = Poi.laptime, aes(x=x, y=y, color = "blue"), size =3) + 
  geom_point(data = Log.laptime, aes(x=x, y=y, color = "red"), size = 3) +
  #geom_point(data = Gam.gvatime, aes(x=x, y=y, color = "green"), size =3) +
  geom_line(data = Poi.laptime, aes(x=x, y=y, color = "blue")) +
  geom_line(data = Log.laptime, aes(x=x, y=y, color = "red")) +
  #geom_line(data = Gam.gvatime, aes(x=x, y=y, color = "green")) +
  scale_color_manual(values =c("blue","red"),labels = c("Poisson","logistic")) +
  scale_y_continuous(
    trans = "log10",                # 纵轴对数变换
    limits = c(5e-2, 5e2),          # 设置纵轴范围
    breaks = c(1e-1, 1e0, 1e1, 1e2, 1e3),  # 设置刻度
    labels = log10_labels           # 使用 10^x 格式标签
  ) +
  scale_x_continuous(
    trans = "log10",                # 横轴对数变换
    limits = c(1e2, 2e5),           # 设置横轴范围
    breaks = c(1e2, 1e3, 1e4, 1e5), # 设置刻度
    labels = log10_labels           # 使用 10^x 格式标签
  ) +
  labs(
    title = "CPU sec. Laplace",     # 图表标题
    x = "Sample size",              # 横轴标签
    y = "Time/s",                   # 纵轴标签
    color = "Method"                # 图例标题
  ) +
  theme_minimal() + 
  theme(legend.position = "top")
LapCT1 <- plot1

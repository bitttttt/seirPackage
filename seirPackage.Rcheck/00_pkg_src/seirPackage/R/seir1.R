library(devtools)
has_devel()
#' ODE solution
#'
#' This function solves an ODE system based on the given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export
model <- function(t, y, param) {
  #设定参数
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- param["R"]
  #设定参数
  beta <- param["beta"]
  mu <- param["mu"]
  gamma <- param["gamma"]
  lamda <- param["lamda"]
  #传染病数学模型
  dSt <- mu * (R - S) - beta * S * I1/R
  dEt <- beta * S * I1/R - mu * E-lamda*E
  dI1t <-   - (mu + gamma) * I1+lamda*E
  dI2t <- gamma * I1 - mu * I2
  #求解结果整合成向量表达
  outcome <- c(dSt, dEt,dI1t, dI2t)
  #返回常微分方程系统求解结果
  list(outcome)
}
#设置评估参数的初始值
times <- seq(0, 156, by = 1/7)
param <- c(mu = 0.000, lamda = 0.04, beta = 5, gamma = 0.2,R = 2)
init <- c(S = 0.9999, E = 0.00009,I1 = 0.00003, I2 = 0)

#调用常微分方程求解函数，传入初始条件， 评估时间，模型以及参数信息
result <- deSolve::ode(y=init, times=times, func=model, parms = param)
result <- as.data.frame(result)

tail(round(result, 3.6),10)


#结果画图
#' @export
seirplot <- ggplot2::ggplot(data=result)+
  ggplot2::geom_line(ggplot2::aes(x=time, y=S,col="S"), lwd=2) +
  ggplot2::geom_line(ggplot2::aes(x=time, y=I1,col="I1"), lwd=2) +
  ggplot2::geom_line(ggplot2::aes(x=time, y=I2,col="I2"), lwd=2) +
  ggplot2::geom_line(ggplot2::aes(x=time, y=E,col="E"), lwd=2) +
  ggplot2::labs(x = "Time",y = "Ratio")+
  ggplot2::scale_color_manual(name = "SEIR",
                              values = c("S" = "pink", "E" = "blue", "I1" ="red", "I2" = "green" ))
#绘制仿真结果并保存为矢量文件
seirplot
ggplot2::ggsave(seirplot, file="seir.pdf", width=7, height=6)
ggplot2::ggsave(seirplot, file="seir.svg", width=7, height=6)

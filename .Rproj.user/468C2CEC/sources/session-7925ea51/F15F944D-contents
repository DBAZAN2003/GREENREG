library(CHAPIREG)
data("datos_anomalia_temperatura")
modelo_MAchapi <- modelo_ma(datos_anomalia_temperatura$anomalia_c, q = 1, include_mean = TRUE)
modelo_MAchapi
plot(modelo_MAchapi)

?modelo_ma


library(stats)
data("datos_anomalia_temperatura")
ts_data <- datos_anomalia_temperatura$anomalia_c
modelo_base <- stats::arima(ts_data,
                            order = c(0, 0, 1),    # p=0, d=0, q=1
                            method = "CSS",
                            include.mean = TRUE)
print(modelo_base)
cat("\nLjung-Box (Residuos):\n")
print(Box.test(residuals(modelo_base), lag = 10, type = "Ljung-Box"))
cat("\nShapiro-Wilk (Residuos):\n")
print(shapiro.test(residuals(modelo_base)))

library(CHAPIREG)
data("datos_nivel_presa")
ts_data <- datos_nivel_presa$nivel_m
modelo_ARIMAchapi <- modelo_arima(ts_data, p = 1, d = 1, q = 1)
modelo_ARIMAchapi
plot(modelo_ARIMAchapi)

?modelo_arima

library(stats)
data("datos_nivel_presa")
ts_data <- datos_nivel_presa$nivel_m
modelo_base <- stats::arima(ts_data,
                            order = c(1, 1, 1),
                            method = "CSS",
                            include.mean = TRUE)
modelo_base
cat("\nLjung-Box (Residuos):\n")
print(Box.test(residuals(modelo_base), lag = 10, type = "Ljung-Box"))
cat("\nShapiro-Wilk (Residuos):\n")
print(shapiro.test(residuals(modelo_base)))

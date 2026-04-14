library(CHAPIREG)
data("datos_roya_cafe")
head(datos_roya_cafe)

modelo_log_chap <- logistico(presencia_roya ~ humedad_relativa_pct + altitud_msnm + manejo_sombra, data = datos_roya_cafe)
modelo_log_chap
plot(modelo_log_chap)

?logistico



library(pROC)
library(ResourceSelection)
data("datos_roya_cafe")

modelo_base <- glm(presencia_roya ~ .,
                    data = datos_roya_cafe,
                    family = binomial(link = "logit"))
summary(modelo_base)
print(exp(coef(modelo_base)))

cat("Prueba Hosmer-Lemeshow")
hl_test <- hoslem.test(modelo_base$y, fitted(modelo_base), g = 10)
print(hl_test)

cat("Curva ROC / AUC")
roc_obj <- roc(modelo_base$y, fitted(modelo_base), quiet = TRUE)
print(auc(roc_obj))

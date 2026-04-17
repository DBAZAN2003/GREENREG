library(CHAPIREG)
data("datos_rendimiento_maiz")
?analisis_datos

maiz<- analisis_datos(datos_rendimiento_maiz, "rendimiento_maiz_ton_ha")

roya<- analisis_datos(datos_roya_cafe, "presencia_roya")

talas<- analisis_datos(datos_tala_ilegal_poisson, "numero_talas")

presa<- analisis_datos(datos_nivel_presa, "nivel_m" )

temp<- analisis_datos(datos_anomalia_temperatura, "anomalia_c")


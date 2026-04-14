## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

# --------------------------------------------------
# Generación de datasets del paquete CHAPIREG
# --------------------------------------------------

# 1. Leer CSV actualizados
datos_anomalia_temperatura <- read.csv(
  "data-raw/datos_anomalia_temperatura.csv"
)

datos_nivel_presa <- read.csv(
  "data-raw/datos_nivel_presa.csv"
)

datos_rendimiento_maiz <- read.csv(
  "data-raw/datos_rendimiento_maiz.csv"
)

datos_roya_cafe <- read.csv(
  "data-raw/datos_roya_cafe.csv"
)

datos_tala_ilegal_poisson <- read.csv(
  "data-raw/datos_tala_ilegal_poisson.csv"
)

# 2. Guardar datasets en el paquete
usethis::use_data(
  datos_anomalia_temperatura,
  datos_nivel_presa,
  datos_rendimiento_maiz,
  datos_roya_cafe,
  datos_tala_ilegal_poisson,
  overwrite = TRUE
)

## GLmulti

datos <- read.csv("datos/data_galicia.txt", sep = " ")


datos$InfFasc <- as.factor(datos$InfFasc)
datos$Aptitude <- as.factor(datos$Aptitude)


glmulti.out <- glmulti(InfFasc ~ Age + Aptitude + Temperature + Rainfall + Altitude + Slope + Density, 
                       family = binomial("logit"), 
                       crit = "aic",
                       data = datos)

saveRDS("./objetos_r/glmulti.rds")
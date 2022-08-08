library(GGally)
library(dplyr)
library(tidyr)
library(mgcv)
library(psych)
library(knitr)
library(R2WinBUGS)
library(caret)
library(pROC)
library(ResourceSelection)

# Carga de los datos
datos <- read.csv("data_galicia.txt", sep = " ")
datos$Aptitude <- as.factor(datos$Aptitude)
datos$InfFasc <- as.factor(datos$InfFasc)

attach(datos)
summary(datos[, c(1, 2, 4, 6, 7, 8, 9, 10)])
summary(datos[, -c(1, 2, 4, 6, 7, 8, 9, 10)])


ggpairs(datos)

par(mfrow=c(3,2))
plot(as.numeric(Age), jitter(as.numeric(InfFasc)-1, factor = 0.2),
     xlab = "Edad",
     ylab = "Presencia del parásito (0 = No, 1 = Sí)")
tt <- table(InfFasc, Age)
prc <- tt[2,] / (tt[1,] + tt[2,])
lines(as.numeric(colnames(tt)), prc,type='b', col='RED')

plot(as.numeric(Temperature), jitter(as.numeric(InfFasc)-1, factor = 0.2),
     xlab = "Temperatura",
     ylab = "Presencia del parásito (0 = No, 1 = Sí)")
tt <- table(InfFasc, Temperature)
prc <- tt[2,] / (tt[1,] + tt[2,])
lines(as.numeric(colnames(tt)), prc,type='b', col='RED')

plot(as.numeric(Rainfall), jitter(as.numeric(InfFasc)-1, factor = 0.2),
     xlab = "Pluviometría",
     ylab = "Presencia del parásito (0 = No, 1 = Sí)")
tt <- table(InfFasc, Rainfall)
prc <- tt[2,] / (tt[1,] + tt[2,])
lines(as.numeric(colnames(tt)), prc,type='b', col='RED')

plot(as.numeric(Altitude), jitter(as.numeric(InfFasc)-1, factor = 0.2),
     xlab = "Altitud",
     ylab = "Presencia del parásito (0 = No, 1 = Sí)")
tt <- table(InfFasc, Altitude)
prc <- tt[2,] / (tt[1,] + tt[2,])
lines(as.numeric(colnames(tt)), prc,type='b', col='RED')

plot(as.numeric(Slope), jitter(as.numeric(InfFasc)-1, factor = 0.2), 
     xlab = "Pendiente",
     ylab = "Presencia del parásito (0 = No, 1 = Sí)")
tt <- table(InfFasc, Slope)
prc <- tt[2,] / (tt[1,] + tt[2,])
lines(as.numeric(colnames(tt)), prc,type='b', col='RED')

plot(as.numeric(Density), jitter(as.numeric(InfFasc)-1, factor = 0.2),
     xlab = "Densidad",
     ylab = "Presencia del parásito (0 = No, 1 = Sí)")
tt <- table(InfFasc, Density)
prc <- tt[2,] / (tt[1,] + tt[2,])
lines(as.numeric(colnames(tt)), prc,type='b', col='RED')

ajuste.logit <- glm(InfFasc ~ ., data=datos, family = binomial(link="logit"))
ajuste.probit <- glm(InfFasc ~ ., data=datos, family = binomial(link="probit"))
ajuste.cloglog <- glm(InfFasc ~ ., data=datos, family = binomial(link="cloglog"))

x <- matrix(c(AIC(ajuste.logit), deviance(ajuste.logit),
              AIC(ajuste.probit), deviance(ajuste.probit),
              AIC(ajuste.cloglog), deviance(ajuste.cloglog)), nrow=3,byrow=T)
colnames(x) <- c("AIC","Deviance")
rownames(x) <- c("ajuste.logit", "ajuste.probit", "ajuste.cloglog")
x

hoslem.test(as.numeric(as.character(datos$InfFasc)), fitted(ajuste.logit))

maximo <- glm(InfFasc ~ ., family = binomial("logit"), data = datos)
minimo <- glm(InfFasc ~ 1, family = binomial("logit"), data = datos)
ajuste.step1 <- step(maximo, direction = "backward", scope = minimo)
ajuste.step1$coefficients

ajuste.step <- glm(InfFasc ~ Y + Aptitude + Temperature + Altitude, data=datos,
            family = binomial(link="logit"))
summary(ajuste.step)


modelo1 <- function(){
  for(i in 1:N){
  InfFasc[i] ~ dbern(p[i])
    logit(p[i]) <- beta1 + beta2*Temperature[i] + beta.A[Aptitude[i]]
  }
  # distribuciones previas
  beta1 ~ dflat()
  beta2 ~ dflat()
  beta.A[2] ~ dflat()
  #restricción
  beta.A[1] <- 0
}
set.seed(1)
datos1 <- list(InfFasc=datos$InfFasc, Temperature=datos$Temperature,
              Aptitude=as.numeric(datos$Aptitude),
              N=dim(datos)[1])
iniciales1 <- function(){
  list(beta1=rnorm(1,0,3), beta2=rnorm(1), beta.A=c(NA,rnorm(1)))}
parametros1 <- c("beta1", "beta2", "beta.A")
ResulModelo1 <- bugs(model = modelo1, data = datos1, inits = iniciales1,
                     param = parametros1, n.iter = 20000, n.burnin = 2000)



ajuste.gam1 <- gam(InfFasc ~ s(X,Y) + Aptitude + te(Altitude, Rainfall), data=datos, 
                  family = binomial(link="logit"))
summary(ajuste.gam1)

par(mfrow = c(2,2))
gam.check(ajuste.gam1)

modelo <- glm(InfFasc ~ Aptitude + Temperature,
              data = datos,
              family=binomial(link="logit"))
summary(modelo)

# Análisis de la capacidad predictiva

# Dividimos los datos
set.seed(123)
split <- sample(c(rep(0, 0.7 * nrow(datos)), rep(1, 0.3 * nrow(datos))))
train <- datos[split == 0, ] 
test <- datos[split == 1, ]

# Estadísticos de predicción para el GLM
train.glm <- glm(InfFasc ~ Y + Aptitude + Temperature + Altitude, data=train, 
                 family = binomial(link="logit"))

p.glm <- predict(train.glm, newdata=test, type='response')
result.glm <- as.factor(ifelse(p.glm > 0.5, 1, 0))
confusionMatrix(data = result.glm, reference = test$InfFasc)
curv_roc.glm <- roc(as.factor(test$InfFasc), p.glm)

# Estadísticos de predicción para el GAM
train.gam <- gam(InfFasc ~ s(X,Y) + Aptitude + te(Altitude, Rainfall), data=train, 
                 family = binomial(link="logit"))

p.gam <- predict(train.gam, newdata=test, type='response')
result.gam <- as.factor(ifelse(p.gam > 0.5, 1, 0))
confusionMatrix(data=result.gam, reference = test$InfFasc)
curv_roc.gam <- roc(as.factor(test$InfFasc), p.gam)

# Representación curvas ROC
par(mfrow=c(1, 2))
plot(curv_roc.glm, xlim=c(1,0), main = "GLM")
plot(curv_roc.gam, xlim=c(1,0), main = "GAM")

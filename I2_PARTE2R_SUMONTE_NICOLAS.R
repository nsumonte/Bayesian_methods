#Parte D Interrogación 2.
library(distributions3)
n_1 = 14
n_2 = 5
y_1prom = 12.56
s1 = 24.65
y_2prom = 17.31
s2 = 11.01
mu_10 = 10
mu_20 = 10
k_10 = 1
k_20 = 5
v_0 = 5
sigma_0 = 10
m1_n1 = (mu_10*k_10 + n_1*y_1prom)/(n_1+k_10)
m2_n2 = (mu_20*k_20 + n_2*y_2prom)/(n_2+k_20)
k1_n1 = n_1 + k_10
k2_n2 = n_2 + k_20
v_n = v_0 + n_1 + n_2
sigma_n = (v_0*sigma_0 + (n_1-1)*s1 + (n_2-1)*s2 + (((k_10*n_1 /(k_10 + n_1)))*(y_1prom-mu_10)^2) + (((k_20*n_2 /(k_20 + n_2)))*(y_2prom-mu_20)^2))/(v_n) 
sigma_n
#D.1
#Hacemos la distribucion de student con vn grados de libertad
T_student_1 <- StudentsT(df = v_n)

#calculamos cota inferior
(m1_n1-m2_n2) + quantile(T_student_1, 0.05 / 2) * (sigma_n) / sqrt(k1_n1+k2_n2)

#calculamos cota superior
(m1_n1-m2_n2) + quantile(T_student_1,1 - 0.05 / 2) * (sigma_n) / sqrt(k1_n1+k2_n2)

#Ahora hacemos lo propio para la priori,  sin calcular los datos, nos queda
T_student_2 <-StudentsT(df = v_0)
#calculamos la cota inferior
(mu_10 - mu_20) + quantile(T_student_2,0.05 / 2) * (2*sigma_0) / sqrt(k_10 + k_20)

#calculamos la cota superior
(mu_10 - mu_20) + quantile(T_student_2,1 - 0.05 / 2) * (2*sigma_0) / sqrt(k_10 + k_20)
#notemos que esta ultima es simetrica en torno al cero y la posteriori no lo es.

#D.2
T_student_3 <- StudentsT(df = n_1 + n_2 - 2)

#calculamos cota inferior
(y_1prom-y_2prom) + quantile(T_student_3, 0.05 / 2) * (s1+s2) / sqrt(n_1+n_2)

#calculamos cota superior
(y_1prom-y_2prom) + quantile(T_student_3,1 - 0.05 / 2) * (s1+s2) / sqrt(n_1+n_2)

#En este caso no se pueden calcular el intervalo de credibilidad para la priori, debido a que v0 tiende a  cero y la chi cuadrado inversa no estaria bien definida

#D.3
#Notemos que las muestras si son muy distintas, esto se debe a que la diferencia de cuadrados es muy distinta y es preponderante al momento de calular los intervalos, además de que el tamaño de diferencia
#es muy distinto en ambas distribuciones.


#PARTE E
#USARÉ MIL MUESTRAS PARA LA SIMULACIÓN
M <- 1000
sigma.p <- 1/rgamma(M,v_n/2,(v_n*sigma_n)/2)
mu_nuevos <- matrix(0,ncol=2,nrow=M)
set.seed(2)
for(i in 1:M){
  mu_nuevos[i,1] <- rnorm(1,m1_n1,sigma.p[i]/k1_n1)
  mu_nuevos[i,2] <- rnorm(1,m2_n2,sigma.p[i]/k2_n2)	
}
#Ahora hacemos juntamos la matriz de mu con el valor de sigma para obtener la tripleta de variables desconocidas.
tres_estimadores = cbind(mu_nuevos,sigma.p)
#Ahora para obtener las y rep, nos basta tomar los datos de las primeras n1 y n2 muestras de nuestra matriz
y_rep1 = matrix(0,ncol=1,nrow=n_1)
y_rep2 = matrix(0,ncol=1,nrow=n_2)
set.seed(2)
for(i in 1:n_1){
  y_rep1[i,1] <- rnorm(1,tres_estimadores[i,1],tres_estimadores[i,3])
}
for(i in 1:n_2){
  y_rep2[i,1] <- rnorm(1, tres_estimadores[i,2],tres_estimadores[i,3])
}
#LUEGO y_rep1 y y_rep2 representan los datos replicados con numero de muestra igual al obtenido.
y_rep3 = matrix(0,ncol=1,nrow=M)
y_rep4 = matrix(0,ncol=1,nrow=M)
for(i in 1:M){
  y_rep3[i,1] <- rnorm(1,tres_estimadores[i,1],tres_estimadores[i,3])
  y_rep4[i,1] <- rnorm(1, tres_estimadores[i,2],tres_estimadores[i,3])
}
#PARTE F
#simulamos las mil muestras
y_rep3 = matrix(0,ncol=2,nrow=M)
for(i in 1:M){
  y_rep3[i,1] <- rnorm(1,tres_estimadores[i,1],tres_estimadores[i,3])
  y_rep3[i,2] <- rnorm(1, tres_estimadores[i,2],tres_estimadores[i,3])
}
#definimos las medidas de discrepancia
Ty1rep <- function(y,...){
  s1 <- var(y[,1])
  s2 <- var(y[,2])
  aux <- log(s1/s2)
  return(aux)
}

Ty2rep <- function(y,M...){
  media1 <- mean(y[,1])
  media2 <- mean(y[,2])
  suma1 <- (sum(y[,1]-media1)^2)/M
  suma2 <- (sum(y[,2]-media2)^2)/M
  aux <- log(suma1/suma2)
  return(aux)
}
Ty1_rep <- numeric()
Ty2_rep <- numeric()

#y notemos que las medias de discrepancia de los datos observados son
T1obs = log((s1/(n_1-1))/(s2/n_2-1))
T2obs = log(s1/s2) 

for(i in 1:M){
  Ty1_rep[i] <- Ty1rep(y_rep3)
  Ty2_rep[i] <- Ty2rep(y_rep3)
}
#finalmente calculamos el valor p predictivo para ambos estimadores

mean(T1obs>Ty1_rep)
mean(T2obs>Ty2_rep)

#si es logico igualar varianzas, esto es debido a que ambos valores p  son muy superiores.




# Chargement des packages nécessaires
library(readxl)
library(lubridate)
library(openxlsx2)
# for the Guinea
# first question
library(deSolve)

# Definition of the function equations differentiates model SEIR
seir_model = function(t, y, params) {
  with(as.list(c(y,params )), {
    dS = -beta * S * I/ N
    dE = beta * S * I / N - sigma * E
    dI = sigma * E - gamma * I
    dR =(1-mu)* gamma * I
    return(list(c(dS, dE, dI, dR)))
  })
}
# C. L. Althaus. Estimating the reproduction number of ebola virus (ebov) during the 2014 outbreak in west
#africa. PLoS currents, 6, 2014
#Données initiales 
N=10^6 # population totale
y_initiales=c(S=N-1, E=0,I=1,R=0)
param=c(beta=1.51*1/5.61,sigma=1/5.3,gamma=1/5.61,mu=0.74)# beta=R0/gamma, Ro=1.51
param
# Temps de simulation

t=seq(0,200,by=1)
# Resolution numerique des equations diff
resol=ode(y=y_initiales,times=t,func = seir_model,parms=param)
resol

# Extraction des valeurs des compartiments S, E, I et R depuis le résultat
S <- resol[, "S"]
E <- resol[, "E"]
I <- resol[, "I"]
R <- resol[, "R"]
par(mfrow=c(2,2))
plot(t, S, type = "l", col = "blue", xlab = "Temps", ylab = "Nombre de personnes", main = "S")
plot(t, E, type = "l", col = "red", xlab = "Temps", ylab = "Nombre de personnes", main = "E")
plot(t, I, type = "l", col = "green", xlab = "Temps", ylab = "Nombre de personnes", main = "I")
plot(t, R, type = "l", col = "black", xlab = "Temps", ylab = "Nombre de personnes", main = "R")

# 3 Tester la sensibilité au taux de transmission
beta_values = c( 0.5, 0.95)
par(mfrow=c(3,4))
for (beta_value in beta_values) {
  resol <- ode(y = y_initiales, times = t, func = seir_model, parms =c(beta=beta_value,sigma=1/5.3,gamma=1/5.61,mu=0.74))
  #resol <- out[, 3]
  S <- resol[, "S"]
  E <- resol[, "E"]
  I <- resol[, "I"]
  R <- resol[, "R"]
  
  # Tracer les courbes
  plot(t, S, type = "l", col = "blue", xlab = "Temps", ylab = "Nombre de personnes", main = "S")
  plot(t, E, type = "l", col = "red", xlab = "Temps", ylab = "Nombre de personnes", main = "E")
  plot(t, I, type = "l", col = "green", xlab = "Temps", ylab = "Nombre de personnes", main = "I")
  plot(t, R, type = "l", col = "black", xlab = "Temps", ylab = "Nombre de personnes", main = "R")
  # Légende
  #legend("topright", c("beta = 0.25", "beta = 0.5", "beta = 0.75"), col = c(rgb(0, 0, 1, 0.25), rgb(0, 0, 1, 0.5), rgb(0, 0, 1, 0.75)), lwd = 2)
  
  
}


# Question 2


library(ggplot2)

# Extraction des données
data = read_xlsx("previous-case-counts.xlsx", sheet = 1)  
data = data[1:265,]
data[,1] = as.Date(data[,1], '%Y-%m-%d')
data_1 = data[data[,1] < as.Date("2015-11-01"),]
Gui_data = data_1[,c(1,2,3)]
Gui_data = data.frame("day"=Gui_data[,1], "TC"=Gui_data[,2], "TD"=Gui_data[,3])

# Tracer TC et TD en fonction du temps
ggplot(data = Gui_data, aes(x = day)) +
  geom_line(aes(y = TC, color = "Total Cases")) +
  geom_line(aes(y = TD, color = "Total Deaths")) +
  labs(x = "Date", y = "Count", color = "Variables") +
  theme_minimal() +
  ggtitle("Évolution des cas et des décès au fil du temps")

# Extraction des données observées
D_obs_TC <- Gui_data$TC
D_obs_TD <- Gui_data$TD
D_obs_TD
# Définition de la fonction de distance
distance <- function(D_obs_TC, D_obs_TD, D_sim_TC, D_sim_TD) {
  dist_TC <- sqrt(sum((D_obs_TC - D_sim_TC)^2))
  dist_TD <- sqrt(sum((D_obs_TD - D_sim_TD)^2))
  return(list(dist_TC = dist_TC, dist_TD = dist_TD))
}


# Définition de la fonction de simulation du modèle SEIR
simulate_SEIR <- function(beta, sigma, gamma, mu, N, T) {
  # Initialisation des compartiments
  S <- N - 1
  E <- 0
  I <- 1
  R <- 0
  # Vecteurs pour stocker les données simulées
  cas_confirmes <- numeric(T)
  deces <- numeric(T)
  # Simulation du modèle SEIR
  for (t in 1:T) {
    # Calcul des flux entre les compartiments
    dS <- -beta * S * I / N
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- (1 - mu) * gamma * I
    # Mise à jour des compartiments
    S <- S + dS
    E <- E + dE
    I <- I + dI
    R <- R + dR
    # Stockage des données simulées
    cas_confirmes[t] <- sigma * E
    deces[t] <- mu * gamma * I
  }
  return(list(cas_confirmes = cas_confirmes, deces = deces))
}


# Algorithme ABC
n <- 1000 # Taille de l'échantillon souhaitée
parametres_acceptes <- list()
# Initialisation des compteurs pour chaque paramètre
compteurs = list(beta = 0, sigma = 0, gamma = 0, mu = 0)

for (i in 1:n) {
  # Simulation des paramètres
  beta <- runif(1, 0, 1)
  sigma <- runif(1, 0, 1)
  gamma <- runif(1, 0, 1)
  mu <- runif(1, 0, 1)
  
  # Simulation des données SEIR
  N <- 10^6  # Population totale
  T <- length(D_obs_TC)
  D_i <- simulate_SEIR(beta, sigma, gamma, mu, N, T)
  D_i_TC <- D_i$cas_confirmes
  D_i_TD <- D_i$deces
  
  # Calcul de la distance entre les données observées et simulées
  dist <- distance(D_obs_TC, D_obs_TD, D_i_TC, D_i_TD)
  dist
  # Seuil de distance
  epsilon <- 50000
  
  # Acceptation des paramètres si la distance est inférieure au seuil
  if (!is.na(dist$dist_TC) & !is.na(dist$dist_TD) &
      dist$dist_TC < epsilon & dist$dist_TD < epsilon) {
    parametres_acceptes[[i]] <- list(beta = beta, sigma = sigma, gamma = gamma, mu = mu)
    # Mise à jour des compteurs pour chaque parametre
    compteurs$beta = compteurs$beta + 1
    compteurs$sigma = compteurs$sigma + 1
    compteurs$gamma = compteurs$gamma + 1
    compteurs$mu = compteurs$mu + 1
  }
}

# Affichage des paramètres acceptés
print(parametres_acceptes)
# Affichage du nombre de chaque parametre accepte
print(compteurs)
  
# Extraction des paramètres acceptés dans un dataframe
param_df =do.call(rbind, parametres_acceptes)
param_df
# Affichage des histogrammes des distributions pour chaque paramètre
par(mfrow=c(2, 2))  # Pour afficher les graphiques dans une grille 2x2
param_df_beta=param_df[,1]
param_df_sigma=param_df[,2]
param_df_gamma=param_df[,3]
param_df_mu=param_df[,4]
param_df_beta
param_df_beta = as.numeric(param_df_beta)
param_df_sigma = as.numeric(param_df_sigma)
param_df_gamma = as.numeric(param_df_gamma)
param_df_mu = as.numeric(param_df_mu)


hist(param_df_beta, main="Distribution de beta", xlab="Valeur de beta")
hist(param_df_sigma, main="Distribution de sigma", xlab="Valeur de sigma")
hist(param_df_gamma, main="Distribution de gamma", xlab="Valeur de gamma")
hist(param_df_mu, main="Distribution de mu", xlab="Valeur de mu")
# Vérifier le type de données des paramètres






# Trier les valeurs des paramètres
sorted_beta <- sort(param_df_beta)
sorted_sigma <- sort(param_df_sigma)
sorted_gamma <- sort(param_df_gamma)
sorted_mu <- sort(param_df_mu)

# Calculer les quantiles à 2.5% et 97.5%
quantiles_beta <- quantile(sorted_beta, c(0.025, 0.975))
quantiles_beta
quantiles_sigma <- quantile(sorted_sigma, c(0.025, 0.975))
quantiles_gamma <- quantile(sorted_gamma, c(0.025, 0.975))
quantiles_mu <- quantile(sorted_mu, c(0.025, 0.975))
# Extraire les bornes des intervalles de crédibilité
 beta_a<- quantiles_beta[1]
beta_b <- quantiles_beta[2]
sigma_a<- quantiles_sigma[1]
sigma_b <- quantiles_sigma[2]
gamma_a<- quantiles_gamma[1]
gamma_b <- quantiles_gamma[2]
mu_a<- quantiles_mu[1]
mu_b <- quantiles_mu[2]

# Afficher les intervalles de crédibilité
cat("Intervalle de crédibilité pour beta : [", beta_a, ", ", beta_b, "]\n")
cat("Intervalle de crédibilité pour sigma : [", sigma_a, ", ", sigma_b, "]\n")
cat("Intervalle de crédibilité pour gamma : [", gamma_a, ", ", sigma_b, "]\n")
cat("Intervalle de crédibilité pour mu : [", mu_a, ", ", mu_b, "]\n")

#Dans cette version mise à jour de la fonction seir_model, nous introduisons la variation de beta(t)
#dans le temps en fonction du paramètre tau qui représente le moment où les mesures de contrôle 
#sont mises en place, et du paramètre k qui représente le taux de décroissance de beta(t)


#question 1 modele amélioré
para=c(beta=1.51*1/5.61,sigma=1/5.3,gamma=1/5.61,mu=0.74,tau=2015,k=0.99)
seir_model_ameliore <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Calcul de beta(t) en fonction du temps t
    beta_t <- ifelse(t <= tau, beta, beta * exp(-k * (t - tau)))
    
    dS <- -beta_t * S * I / N
    dE <- beta_t * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- (1 - mu) * gamma * I
    
    return(list(c(dS, dE, dI, dR)))
  })
}
resol_ameliore=ode(y=y_initiales,times=t,func = seir_model_ameliore,parms=para)
resol_ameliore

# Extraction des valeurs des compartiments S, E, I et R depuis le résultat
S <- resol_ameliore[, "S"]
E <- resol_ameliore[, "E"]
I <- resol_ameliore[, "I"]
R <- resol_ameliore[, "R"]
par(mfrow=c(2,2))
plot(t, S, type = "l", col = "blue", xlab = "Temps", ylab = "Nombre de personnes", main = "S")
plot(t, E, type = "l", col = "red", xlab = "Temps", ylab = "Nombre de personnes", main = "E")
plot(t, I, type = "l", col = "green", xlab = "Temps", ylab = "Nombre de personnes", main = "I")
plot(t, R, type = "l", col = "black", xlab = "Temps", ylab = "Nombre de personnes", main = "R")


#question 3 deuxième partie
#Du debut de l'epidemie jusqu'au 13/04/2016
data=data[1:265,]
data
data[,1] = as.Date(data[,1], '%Y-%m-%d')
data_2 = data[data[,1] < as.Date("2016-04-14"),]
data_2
Gui_data = data_2[,c(1,2,3)]
Gui_data=data.frame("day"=Gui_data[,1],"TC"=Gui_data[,2],"TD"=Gui_data[,3])
Gui_data
library(ggplot2)

# Convertir la colonne 'day' en un objet de date
Gui_data$day =as.Date(Gui_data$day)

# Tracer TC et TD en fonction du temps
ggplot(data = Gui_data, aes(x = day)) +
  geom_line(aes(y = TC, color = "Total Cases")) +
  geom_line(aes(y = TD, color = "Total Deaths")) +
  labs(x = "Date", y = "Count", color = "Variables") +
  theme_minimal() +
  ggtitle("Évolution des cas et des décès au fil du temps")




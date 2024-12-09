# Inclure le fichier qui contient le module MarketSimulation
include("MarketSimulation.jl")

# Charger le module MarketSimulation
using .MarketSimulation  # Utiliser le "." pour indiquer un module local

nb_agents = 1000

# Initialisation de 2 agents
#T, gamma, lambda, a, b, tmin, tmax, pmin, pmax, n_producteur = MarketSimulation.initialize_test()
#println("Agents initialisés avec succès")

T, gamma, lambda, a, b, tmin, tmax, pmin, pmax, n_producteur, n_consommateur, n_mixte = MarketSimulation.initialize_agents(nb_agents)
println("Agents initialisés avec succès")

rho = 10.0
max_iters = 10000

# Exécution de la simulation pour 100 itérations avec un facteur de pénalité ρ = 0.1
MarketSimulation.simulate_market(T, n_producteur, max_iters, a, b, tmin, tmax, pmin, pmax, gamma, lambda, rho/10)

println("Nombre de producteurs : ", n_producteur)
println("Nombre de consommateurs : ", n_consommateur)
println("Nombre de prosumers : ", n_mixte)

println("Simulation terminée")

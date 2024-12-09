using Pkg
Pkg.add("JuMP")
Pkg.add("Plots")
Pkg.add("PGFPlotsX")

using JuMP
using Plots  # Pour la visualisation
using PGFPlotsX  # Backend pour graphiques vectoriels
pgfplotsx()  # Active le backend PGFPlotsX
 
 
module MarketSimulation
    using JuMP
    using Random
    using Plots  # Pour la visualisation
    gr()  # Active le backend GR pour l'affichage
    include("MarketAgents.jl")
    using .MarketAgents: Agent, admm_update, update_lagrange, cost_function

    function simulate_market(
        T::Matrix{Float64}, n_producteur::Int64, max_iters::Int64, 
        a::Vector{Float64}, b::Vector{Float64}, tmin::Vector{Float64}, 
        tmax::Vector{Float64}, pmin::Vector{Float64}, pmax::Vector{Float64}, 
        gamma::Matrix{Float64}, lambda::Matrix{Float64}, rho::Float64
    )
        n, m = size(T)
        R = zeros(n)
        S = zeros(n)
        P = zeros(n)
        Mu = zeros(n)
        T_mean = zeros(n)
        errors = Float64[]  # Liste des erreurs au cours des itérations
        times = Float64[]   # Liste des temps de chaque itération
        max_error = 1e-6  # Tolérance de l'erreur
        error = 2 * max_error  # Initialisation de l'erreur pour débuter la boucle
        k = 0
        bt1 = zeros(n, m)  # Matrice de mise à jour des prix intermédiaires
        n_mixte = 0 # nombre de prosumers à changer

        

        while k < max_iters && error > max_error
            println("Itération $k")
    
            k += 1
            T_old = deepcopy(T)
            # Mesurer le temps de l'itération
            start_time = time()  # Début de l'itération

            # Mise à jour pour chaque agent basé sur son ID correspondant à l'itération
            T_bis, P, Mu, T_mean = admm_update(T, n_producteur, n_mixte, max_iters, a, b, tmin, tmax, pmin, pmax, gamma, lambda, rho, bt1, max_error, P, Mu, T_mean)
            
            # Mise à jour des multiplicateurs de Lagrange
            R, S, lambda, bt1 = update_lagrange(T_old, T_bis, lambda, rho)
            
            #println("T_new: ", T_bis)
            
        
            # Affichage des prix des agents après chaque itération
            #println("Prix des agents à l'itération $k : ", [lambda[i, :] for i in 1:size(lambda, 1)])
            
            println("error S :", maximum(S))
            println("error R :", maximum(R))
            # Calcul de l'erreur maximale entre voisins
            error = max(maximum(S), maximum(R))
            push!(errors, error)
            println("Erreur maximale de prix à l'itération $k: $error")
            println("***************************************************")

             # Mesurer le temps de l'itération
             end_time = time()  # Fin de l'itération
             push!(times, end_time - start_time)  # Ajouter le temps de l'itération à la liste des temps
             

            T = deepcopy(T_bis)
        end
        
         # Tracé de l'évolution de l'erreur en fonction du temps
         iterations = 2:length(errors)  # Abscisses correspondant aux itérations
         plt1 = scatter(iterations, errors[2:end], title="Évolution de l'erreur au cours des itérations", xlabel="Itérations", ylabel="Erreur", legend=false, marker=:circle, color=:blue, xticks=100:100:length(errors), yscale=:log10, xtickfont=(:normal, 10))
         display(plt1)  # Affichage explicite du graphique de l'erreur
         savefig("error_plot_1000_v2.svg")  # Sauvegarde en SVG
         println("Le graphique de l'erreur a été sauvegardé sous 'error_plot.svg'.")
     
         # Calcul manuel de la moyenne
        total_time = sum(times)
        num_iterations = length(times)
        average_time = total_time / num_iterations

        # Tracé du temps de calcul en fonction des itérations
        plt2 = scatter(iterations,times,title="Temps de calcul par itération",xlabel="Itérations",ylabel="Temps (secondes)",legend=false,marker=:circle,color=:red,xticks=100:100:length(times),xtickfont=(:normal, 10))

        # Ajout de la ligne horizontale représentant la moyenne
        hline!([average_time], color=:blue, linestyle=:dash, label="Moyenne : $(round(average_time, digits=3)) secondes")

        # Par exemple, décaler l'annotation verticalement en ajoutant ou soustrayant une valeur
        annotate!(length(iterations) / 2, average_time + 0.005, text("Moyenne : $(round(average_time, digits=3)) s", 10, :blue))


        # Affichage du graphique
        display(plt2)

        # Sauvegarde du graphique
        savefig("time_plot_1000_v2.svg")  # Sauvegarde en SVG
        println("Le graphique du temps de calcul a été sauvegardé sous 'time_plot.svg'.")
        
    end
    
    function simulate_markets(
        T::Matrix{Float64}, n_producteur::Int64, max_iters::Int64, 
        a::Vector{Float64}, b::Vector{Float64}, tmin::Vector{Float64}, 
        tmax::Vector{Float64}, pmin::Vector{Float64}, pmax::Vector{Float64}, 
        gamma::Matrix{Float64}, lambda::Matrix{Float64}, rho::Float64
    )
        n, m = size(T)
        R = zeros(n)
        S = zeros(n)
        P = zeros(n)
        Mu = zeros(n)
        T_mean = zeros(n)
        errors = Float64[]  # Liste des erreurs au cours des itérations
        times = Float64[]   # Liste des temps de chaque itération
        max_error = 1e-6  # Tolérance de l'erreur
        error = 2 * max_error  # Initialisation de l'erreur pour débuter la boucle
        k = 0
        bt1 = zeros(n, m)  # Matrice de mise à jour des prix intermédiaires
        n_mixte = 0 # nombre de prosumers à changer

        

        while k < max_iters && error > max_error
            println("Itération $k")
    
            k += 1
            T_old = deepcopy(T)
            # Mesurer le temps de l'itération
            start_time = time()  # Début de l'itération

            # Mise à jour pour chaque agent basé sur son ID correspondant à l'itération
            T_bis, P, Mu, T_mean = admm_update(T, n_producteur, n_mixte, max_iters, a, b, tmin, tmax, pmin, pmax, gamma, lambda, rho, bt1, max_error, P, Mu, T_mean)
            
            # Mise à jour des multiplicateurs de Lagrange
            R, S, lambda, bt1 = update_lagrange(T_old, T_bis, lambda, rho)
    

            #println("T: ", T_old)
            #println("T: ", T_bis)
            
        
            # Affichage des prix des agents après chaque itération
            #println("Prix des agents à l'itération $k : ", [lambda[i, :] for i in 1:size(lambda, 1)])
            
            println("error S :", maximum(S))
            println("error R :", maximum(R))
            # Calcul de l'erreur maximale entre voisins
            error = max(maximum(S), maximum(R))
            push!(errors, error)
            println("Erreur maximale de prix à l'itération $k: $error")
            println("***************************************************")

             # Mesurer le temps de l'itération
             end_time = time()  # Fin de l'itération
             push!(times, end_time - start_time)  # Ajouter le temps de l'itération à la liste des temps
             

            T = deepcopy(T_bis)
        end
     
         # Calcul manuel de la moyenne
        total_time = sum(times)
        num_iterations = length(times)
        average_time = total_time / num_iterations
        return average_time
        
    end
    
    
    # Initialisation des agents et du réseau de voisins
    function initialize_agents(n_agents::Int)
        roles = ["producteur", "consommateur", "mixte"]
        lambda = zeros(n_agents, n_agents)  # prix
        gamma = zeros(n_agents, n_agents)  # pénalités
        T = zeros(n_agents, n_agents)  # transferts
        a = zeros(n_agents)  # coefficient pour la demande
        b = zeros(n_agents)  # coefficient pour le prix
        tmax = zeros(n_agents)  # transferts max
        tmin = zeros(n_agents)  # transferts min
        pmax = zeros(n_agents)  # production max
        pmin = zeros(n_agents)  # production min
    
        # Déterminer le nombre d'agents producteurs, consommateurs et mixtes
        n_producteur = div(n_agents, 3)  # Division entière
        #rand(1:n_agents - 2)  # Minimum 1 consommateur et 1 mixte
        n_consommateur = div(n_agents, 3)  # Division entière
        #rand(1:(n_agents - n_producteur - 1))
        n_mixte = n_agents - n_producteur - n_consommateur
    
        for i in 1:n_agents
            if i <= n_producteur
                # Producteur pur
                a[i] = 1
                b[i] = 0.5
                pmax[i] = 15
                pmin[i] = 0
                tmax[i] = pmax[i]
                tmin[i] = 0
                role = roles[1]  # "producteur"
            elseif i <= n_producteur + n_consommateur
                # Consommateur pur
                a[i] = 1
                b[i] = 0.5
                pmax[i] = 0
                pmin[i] = -15
                tmax[i] = 0
                tmin[i] = pmin[i]
                role = roles[2]  # "consommateur"
            else
                # Agent mixte
                a[i] = 1
                b[i] = 0  
                pmax[i] = 15
                pmin[i] = -15
                tmax[i] = pmax[i]
                tmin[i] = pmin[i]
                role = roles[3]  # "mixte"
            end
        end
    
        return T, gamma, lambda, a, b, tmin, tmax, pmin, pmax, n_producteur, n_consommateur, n_mixte
    end
    
      # Initialisation des agents et du réseau de voisins
    function initialize_test()
        n_agents = 2
        n_producteur = 1
        n_mixte = 0
        roles = ["producteur", "consommateur"]
        lambda = zeros(n_agents, n_agents)  # prix
        gamma = zeros(n_agents, n_agents)  # pénalités
        T = zeros(n_agents, n_agents)
        a = zeros(n_agents)
        b = zeros(n_agents)
        tmax = zeros(n_agents)
        tmin = zeros(n_agents)
        pmax = zeros(n_agents)
        pmin = zeros(n_agents)
        
        n_consommateur = n_agents - n_producteur
        
        for i in 1:n_agents
            id = i
            if i<=n_producteur
                a[i] = 1
                b[i] = 4
                pmax[i] = 60
                pmin[i] = 0
                tmax[i] = pmax[i]
                tmin[i] = 0
                role = roles[1]
            elseif i>n_producteur
                a[i] = 1
                b[i] = 8
                pmax[i] = 0
                pmin[i] = -30
                tmax[i] = 0
                tmin[i] = pmin[i]
                role = roles[2]
            end
        end
        return T, gamma, lambda, a, b, tmin, tmax, pmin, pmax, n_producteur, n_consommateur, n_mixte
    end  

end

    
      





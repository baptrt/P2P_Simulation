module MarketAgents

mutable struct Agent
    id::Int64         
end

# Fonction coût locale mise à jour à chaque itération en fonction de P
function cost_function(agent::Agent, a::Vector{Float64}, P::Float64)
    a, b = a[agent.id], b[agent.id]
    return a * P^2 + b * P
end

##Cas pour Producteur-Consommateur
# Mise à jour d'ADMM pour optimiser les trades d'un agent spécifique
#=
function admm_update(
    T::Matrix{Float64}, n_producteur::Int64, max_iters::Int64,
    a::Vector{Float64}, b::Vector{Float64}, tmin::Vector{Float64}, tmax::Vector{Float64}, 
    pmin::Vector{Float64}, pmax::Vector{Float64}, gamma::Matrix{Float64}, 
    lambda::Matrix{Float64}, rho::Float64, bt1::Matrix{Float64}, 
    max_error::Float64, P::Vector{Float64}, Mu::Vector{Float64}, 
    T_mean::Vector{Float64}
)
    # Taille de la matrice T
    n, m = size(T)
    rhol = 1.0
    T_new = T

    for i in 1:n
        iter = 0
        erreur = 2 * max_error
        t_mean, p, mu = 0.0, 0.0, 0.0
        voisins = i <= n_producteur ? (n_producteur + 1) : n : 1:n_producteur
        n_voisin = length(voisins)

        while (erreur > max_error && iter < max_iters) || iter < 2
            iter += 1
            s = 0.0
            list_err = zeros(n)

            for j in voisins
                bt2 = T[i, j] - t_mean + p - mu
                T_new[i, j] = (rho * bt1[i, j] + bt2 * rhol - gamma[i, j]) / (rho + rhol)
                T_new[i, j] = max(min(T_new[i, j], tmax[i]), tmin[i])
                s += T_new[i, j]
                list_err[j] = abs(T_new[i, j] - T[i, j])
                T[i, j] = T_new[i, j]
            end

            erreur = maximum(list_err)
            t_mean = s / n_voisin
            bp1 = mu + t_mean
            p = (n_voisin * rhol * bp1 - n_voisin * b[i]) / (n_voisin * rhol + n_voisin^2 * a[i])
            p = max(min(p, pmax[i] / n_voisin), pmin[i] / n_voisin)
            mu = mu + t_mean - p
            erreur = max(erreur, abs(t_mean - p))
        end

        P[i] = p
        Mu[i] = mu
        T_mean[i] = t_mean
    end
    return T, P, Mu, T_mean
end
=#

##Cas pour Producteur-Consommateur-Mixte
# Mise à jour d'ADMM pour optimiser les trades d'un agent spécifique
function admm_update(
    T::Matrix{Float64}, n_producteur::Int64, n_mixte::Int64, max_iters::Int64,
    a::Vector{Float64}, b::Vector{Float64}, tmin::Vector{Float64}, tmax::Vector{Float64}, 
    pmin::Vector{Float64}, pmax::Vector{Float64}, gamma::Matrix{Float64}, 
    lambda::Matrix{Float64}, rho::Float64, bt1::Matrix{Float64}, 
    max_error::Float64, P::Vector{Float64}, Mu::Vector{Float64}, 
    T_mean::Vector{Float64}
)
    # Taille de la matrice T
    n, m = size(T)
    rhol = 10.0
    T_new = T

    # Indices des groupes
    n_consommateur = n - n_producteur - n_mixte
    prod_indices = 1:n_producteur
    cons_indices = (n_producteur + 1):(n_producteur + n_consommateur)
    mixte_indices = (n_producteur + n_consommateur + 1):n

    for i in 1:n
        iter = 0
        erreur = 2 * max_error
        t_mean, p, mu = 0.0, 0.0, 0.0

        # Définir les voisins selon le rôle de l'agent
        voisins = if i in prod_indices
            cons_indices
        elseif i in cons_indices
            prod_indices
        else  # Agent mixte
            union(prod_indices, cons_indices)
        end

        n_voisin = length(voisins)

        # Mise à jour itérative
        while (erreur > max_error && iter < max_iters) || iter < 2
            iter += 1
            s = 0.0
            list_err = zeros(n)

            for j in voisins
                bt2 = T[i, j] - t_mean + p - mu
                T_new[i, j] = (rho * bt1[i, j] + bt2 * rhol - gamma[i, j]) / (rho + rhol)
                T_new[i, j] = max(min(T_new[i, j], tmax[i]), tmin[i])
                s += T_new[i, j]
                list_err[j] = abs(T_new[i, j] - T[i, j])
                T[i, j] = T_new[i, j]
            end

            erreur = maximum(list_err)
            t_mean = s / n_voisin
            bp1 = mu + t_mean
            p = (n_voisin * rhol * bp1 - n_voisin * b[i]) / (n_voisin * rhol + n_voisin^2 * a[i])
            p = max(min(p, pmax[i] / n_voisin), pmin[i] / n_voisin)
            mu = mu + t_mean - p
            erreur = max(erreur, abs(t_mean - p))
        end

        P[i] = p
        Mu[i] = mu
        T_mean[i] = t_mean
    end

    return T, P, Mu, T_mean
end


# Mise à jour des multiplicateurs de Lagrange
function update_lagrange(T::Matrix{Float64}, T_new::Matrix{Float64}, lambda::Matrix{Float64}, rho::Float64)
    n, m = size(T)
    R = zeros(n)
    S = zeros(n)
    bt1 = zeros(n, m)
    
    for i in 1:n 
        r = 0
        s = 0
        for j in 1:m  # Boucle sur m
            lambda[i, j] = lambda[i, j] + rho / 2 * (T_new[i, j] + T_new[j, i])
            bt1[i, j] = 0.5*(T_new[i,j] - T_new[j,i]) - lambda[i,j]/rho
            r += abs(T_new[i, j] + T_new[j, i])
            s += abs(T_new[i, j] - T[i, j])
        end
        R[i] = r
        S[i] = s
    end
    # Retourner les résultats
    return R, S, lambda, bt1
end

end

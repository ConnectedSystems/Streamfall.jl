NSE(obs, modeled) = 1.0 - sum((obs .- modeled).^2) / sum((obs .- mean(obs)).^2)

NNSE(obs, modeled) = 1.0 / (2.0 - NSE(obs, modeled))

RMSE(obs, modeled) = (sum((modeled .- obs).^2)/length(modeled))^0.5


"""Coefficient of determination (R^2)
"""
function R2(obs::Array, modeled::Array)
    ss_tot = sum((obs .- mean(obs)).^2)
    ss_res = sum((obs .- modeled).^2)

    local r2 = 1 - (ss_res / ss_tot)
    return r2
end


"""Determine adjusted R^2

Parameters
----------
obs : observations
modeled : modeled results
p : number of explanatory variables
"""
function ADJ_R2(obs::Array, modeled::Array, p::Int64)
    n = length(obs)
    adj_r2 = 1 - (1 - R2(obs, modeled)) * ((n - 1) / (n - p - 1))

    return adj_r2
end

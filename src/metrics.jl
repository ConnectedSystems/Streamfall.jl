NSE(obs, modeled) = 1.0 - sum((obs .- modeled).^2) / sum((obs .- mean(obs)).^2)

NNSE(obs, modeled) = 1.0 / (2.0 - NSE(obs, modeled))

RMSE(obs, modeled) = (sum((modeled .- obs).^2)/length(modeled))^0.5

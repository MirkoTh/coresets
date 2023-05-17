B <- 10
d <- 2
k <- 2
sigma_bar_upper <- .5
m_proc <- 10000
s_core <- seq(100, 1000, by = 100)


# epsilon model
plot((B^2*d)/(k^(2/d)))




# epsilon estimation
# upper bound
plot(sigma_bar_upper * B^2 * ((sqrt(k*d))/sqrt(m_proc)))
# lower bound
plot(sigma_bar_upper * B^2 * ((sqrt(k^(1-4/d)))/sqrt(m_proc)))



plot(sqrt(d*k)/(sqrt(s_core) - sqrt(d*k)))

N_vals = c(250, 500, 750)
mu_vals = c(0, 1.5)

sim_frame = expand.grid(N = N_vals, mu = mu_vals)
sim_frame
# 250 0.0
# 500 0.0
# 750 0.0
# 250 1.5
# 500 1.5
# 750 1.5
write.table(sim_frame, file = "inputs.txt", 
            col.names = FALSE, row.names = FALSE)

# Create arguments(parameters) to pass to the job array in Slurm

rm(list = ls())
batch_1 <- c(0, 0, 0)
batch_2 <- c(0.6, 0, 0)
batch_3 <- c(0.7, 0, 0)
batch_4 <- c(0.8, 0, 0)
batch_5 <- c(0.8, 0.1, 0)

# Fixed vaccination column names
vaccination_times <- c("vacc_2", "vacc_4", "vacc_9")

# Dynamically find all batch_* objects in environment
batch_vars <- ls(pattern = "^batch_[0-9]+$")

# Sort batch_vars numerically by batch number
batch_vars <- batch_vars[order(as.numeric(sub("batch_", "", batch_vars)))]

# Retrieve the batch vectors in sorted order
batch_list <- mget(batch_vars)

# Check all batches have the right length
stopifnot(all(sapply(batch_list, length) == length(vaccination_times)))

# Combine into dataframe
df <- as.data.frame(do.call(rbind, batch_list))
rownames(df) <- batch_vars
colnames(df) <- vaccination_times

# Check row sums ≤ 1
if (any(rowSums(df) > 1)) {
  stop("One or more rows have total probability > 1")
}

#print(df)
write.table(df, file = "data/params_job_array/inputs.txt", 
            col.names = FALSE, row.names = FALSE)

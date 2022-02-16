library(knitr)
library(dplyr)
library(ggplot2)
library(ompr)
library(ompr.roi)
library(ROI)
library(ROI.plugin.glpk)

set.seed(54321)

I <- 10
K <- 2

locations <- data.frame(id = 1:I,
                        x=runif(I), 
                        y=runif(I))

travel_costs <- as.matrix(dist(select(locations, x, y), diag = TRUE, upper = TRUE))

depots <- sample(locations$id, K)

job_costs <- rpois(I,1)+1
job_costs[depots] <- 0

locations <- locations %>% mutate(is_depot = id %in% depots)

capacities <- c(15,10)

travel_costs_func <- function(i, i_) {
  vapply(seq_along(i), function(k) travel_costs[i[k], i_[k]], numeric(1L))
}

job_costs_func <- function(i) {
  vapply(seq_along(i), function(k) job_costs[i[k]], numeric(1L))
}

# the depot is always idx 1
model <- MILPModel() %>%
  
  # we create a variable that is 1 if we travel from city i to i_ by Salesman k
  add_variable(x[i, i_, k], 
               i = 1:I,
               i_ = 1:I, 
               k = 1:K, type = "binary") %>%
  
  # minimize travel costs and latest arrival
  set_objective(sum_over(colwise(travel_costs_func(i, i_)) * x[i, i_, k] + 
                           colwise(job_costs_func(i_)) * x[i, i_, k],
                         i = 1:I,
                         i_ = 1:I,
                         k = 1:K), "min") %>%
  
  # Total travel and i_ob costs must not exceed engineer capacity
  add_constraint(sum_over(colwise(travel_costs_func(i, i_)) * x[i, i_, k] +
                            colwise(job_costs_func(i_)) * x[i, i_, k],
                          i = 1:I,
                          i_ = 1:I) <= capacities[k], 
                 k=1:K) %>%
  
  # you cannot go to the same city
  add_constraint(x[i, i, k] == 0, 
                 i = 1:I, 
                 k = 1:K) %>%
  
  # each salesman needs to leave the depot
  add_constraint(sum_over(x[depots[k], i_, k],
                          i_ = 1:I) == 1,
                 k = 1:K) %>%
  
  # each salesman needs to come back to the depot
  add_constraint(sum_over(x[i, depots[k], k], 
                          i = 1:I) == 1,
                 k = 1:K) %>%
  
  # a helper variable for the MTZ formulation of the tsp
  add_variable(u[i, k], i = 1:I, k = 1:K, lb = 1, ub = I) %>%
  
  # if a salesman comes to a city he has to leave it as well
  add_constraint(sum_over(x[i_, i, k], i_ = 1:I) == sum_over(x[i, i_, k], i_ = 1:I), 
               i = 1:I,
               k = 1:K) %>%
  
  # leave each city with only one salesman
  add_constraint(sum_over(x[i, i_, k], i_ = 1:I, k = 1:K) == 1, i = 1:I) %>%
  
  # arrive at each city with only one salesman
  add_constraint(sum_over(x[i, i_, k], i = 1:I, k = 1:K) == 1, i_ = 1:I)


# Add MTZ constraints for each k worker
# ensure no subtours (arc constraints)
for (z in 1:K) {
  model <- model %>%
    add_constraint(u[i, k] >= 2, i = (1:I)[-depots[z]], k = z) %>%
    add_constraint(u[i, k] - u[i_, k] + 1 <= (I - 1) * (1 - x[i, i_, k]), 
                   i = (1:I)[-depots[z]], i_ = (1:I)[-depots[z]], k = z)
}

result <- solve_model(model, with_ROI(solver = "symphony",
                                      verbosity = -2, 
                                      time_limit = 60 * 3))

solution <- get_solution(result, x[i, i_, k]) %>%
  filter(value > 0) %>%
  as_tibble()

solution$travel_cost <- apply(solution, 1, function(r){ travel_costs[as.integer(r["i"]), as.integer(r["i_"])] })
solution$job_cost <- job_costs[solution$i_]

solution%>%
  group_by(k) %>%
  summarise(total_travel = sum(travel_cost),
            total_job = sum(job_cost),
            total_time = total_travel + total_job)

paths <- select(solution, i, i_, k) %>%
  rename(from = i, to = i_) %>%
  mutate(trip_id = row_number()) %>%
  tidyr::gather(property, idx_val, from:to) %>%
  mutate(idx_val = as.integer(idx_val)) %>%
  inner_join(locations, by = c("idx_val" = "id"))

ggplot(locations, aes(x, y)) + 
  geom_point(aes(size = is_depot)) + 
  geom_path(data = paths, aes(group = trip_id, color = factor(k)),
            arrow=arrow(length = unit(0.3, "cm"))) + 
  ggtitle(paste0("Optimal cost: ", round(objective_value(result), 2))) +
  labs(colour= "Engineer", 
       size= "Start Location",
       x=NULL, y=NULL) +
  scale_colour_brewer(type = "qual", palette = 2)


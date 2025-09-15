rm(list = ls())

##load in packages
library(patchwork)
library(performance)
library(DHARMa)
library(fitdistrplus)
library(gamlss)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(mgcv)
library(see)
library(FSA)
library(lme4)
library(car)
library(lmtest)
library(glmmTMB)


#data downloaded from RLS website
#unnecessary columns deleted, non-Bates lab data deleted (KCCA, BMKC, BMSC30-33)
#Swiss Boy BMSC28 (May) changed to Swiss Boy (BMSC24)
#column added for survey_region as Central Coast (Rick surveys 2017), CR, BS, RR, Vic
#add in 2024 data grabbed from google drive final datsheet
#went through all sites in datasheet and identified years where Pycnopodia was absent
## double checked with 2024 dataset and full invert dataset to differentiate between
### site not sampled or pycnopodia absent, year where site was not sampled was
#### ommitted while site & year were added when sampled but pycnopodia absent
#####did this for VIC, BMSC, RR, and CR sites
#number of transects for each site was determined from master downloaded datasheet
#and were then diffirentiated into different transects based on depth
###this was used for abundance -> density measurements

##load data
allPD <- total_pycno_AllPycnoData_220924

pycno_dat <-allPD %>%
  mutate(survey_year = as.factor(survey_year))

View(pycno_dat)

pycno_BSdat <- allPD %>%
  filter(survey_region == "BS") %>%
  mutate(survey_year = as.factor(survey_year))


#seperating forest & understory gives no significance, make into present/absent again
##merge forest and low and understory kelp.pa with absent
pycno_BSdat <- pycno_BSdat %>%
  mutate(kelp.pa = case_when(
    kelp.pa == "low" ~ "present",
    kelp.pa == "understory" ~ "present",
    kelp.pa == "forest" ~ "present",
    TRUE ~ kelp.pa  # Keep other values as they are
  ))


## Q1: Pycnopodia recovery - is there a difference in total between years?
##################### GLMM

#### FIRST USED GLMER TO ACCOUNT FOR SITE AS A RANDOM EFFECT
### QUESTION: Does it make sense to make site a random effect if kelp presence/absence
#is tightly linked to site

boxplot(data=pycno_BSdat, total ~ kelp.pa)

## WORKED FOR TOTAL DIFFS BETWEEN SURVEY YEARS BELOW
mod_pycno1_glmm <- glm(count ~ Transect.level.kelp + as.factor(survey_year) + (1 | site_code),
                         family = poisson,
                         data = pycno_dat)

summary(mod_pycno1_glmm)
#AIC: 444.7 , p2023= 0.001774

#when we add kelp presence/absence as an additional predictor
#NO RELATIONSHIP IS OBSERVED BETWEEN KELP.PA ABSENT AND FOREST/UNDERSTORY OR PRESENT (CLUMPED)
mod_pycno2_glmm <- glmer(total ~ kelp.pa + (1 | site_name),
                         family = poisson,
                         data = pycno_BSdat)

summary(mod_pycno2_glmm)
#AIC: 432.5, p2021= 0.219665, p2022=  0.011656,  p2023 =  0.000195, p.present = 0.537959

#If ignore site as a random effect and use GLM, significance if we clump forest 
#and understory and low into a general present category

# Fit a Poisson model with kelp presence predicting Pycnopodia counts
poisson_model <- glm(total ~ kelp.pa, family = poisson(link = "log"), data = pycno_BSdat)
# Calculate dispersion (residual deviance / degrees of freedom)
dispersion <- sum(residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
dispersion
# data is overdispersed (>1, 3.514441) so fitting it w a Negative Binomial regression
nb_model <- glm.nb(total ~ kelp.pa + depth, data = pycno_BSdat)
summary(nb_model)


# this model gave a warning message re: non-convergence (optimization problem)
#can arise from multiple factors, decided to use he bobyqa optimizer as opposed to the
#default optimizer
mod_pycno3_glmm <- glmer(total ~ comparison_date + kelp.pa + depth + (1 | site_name),
                         family = poisson,
                         data = pycno_BSdat,
                         control = glmerControl(optimizer = "bobyqa"))
summary(mod_pycno3_glmm)
#AIC: 424.9, p2021=  0.49521, p2022= 0.02421, p2023 = 3.88e-05, pdepth = 0.09749
#pkelp.forest= 0.59001 , punderstory= 0.00875, plow=4.85e-05 

bestmodel <- mod_pycno3_glmm


#Check for multicollinearity
library(car)
vif_values <- vif(bestmodel)
print(bestmodel)
#no significant multicollinearity issues among predictors (all under 5)


#WHY POISSON? because count data distributions (e.g., visit counts) often have a 
#Poisson distribution, Poisson regression tends to fit these data better than 
#linear regression does (which assumes a normal distribution) 

#checking homoscedasticity of data in bets fitting model 
#The Breusch-Pagan test checks for heteroscedasticity in the residuals of a 
#linear regression model (variance of residuals/errors is not constant across years)
#H0 = variance of residuals is constant

bptest(bestmodel)
#p-value = 0.1685, fail to reject null, - homoscedastic data!!










#now lets add size_class in
mod_pycno4_glmm <- glmer(total ~ comparison_date + kelp.pa + size_class + depth + (1 | site_name),
                         family = poisson,
                         data = pycno_BSdat,
                         control = glmerControl(optimizer = "bobyqa"))
summary(mod_pycno4_glmm)





mod_pycno4_glmm <- glmer(total ~ comparison_date + kelp.pa + transect + size_class + (1 | site_name),
                         family = poisson,
                         data = pycno_BSdat)
summary(mod_pycno4_glmm)
#AIC: 344.7, BEST model, p2021= 0.35151, p2022= 0.80506, p2023 = 0.01039
#Need to rescale variables to account for the large eigenvalue ratio
# Example of standardizing a variable
pycno_BSdat$comparison_date <- as.numeric(as.character(pycno_BSdat$comparison_date))
pycno_BSdat$transect <- as.numeric(as.character(pycno_BSdat$transect))

# Function for min-max scaling
min_max_scale <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

pycno_BSdat$comparison_date <- min_max_scale(pycno_BSdat$comparison_date)
pycno_BSdat$transect <- min_max_scale(pycno_BSdat$transect)

# Example of log transformation (make sure to add a small constant if any values are zero)
pycno_BSdat$size_class <- log(pycno_BSdat$size_class + 1)

# Centering variables
pycno_BSdat$comparison_date <- pycno_BSdat$comparison_date - mean(pycno_BSdat$comparison_date)
pycno_BSdat$transect <- pycno_BSdat$transect - mean(pycno_BSdat$transect)

# Check structure and convert to factor if necessary
if (!is.factor(pycno_BSdat$site_name)) {
  pycno_BSdat$site_name <- as.factor(pycno_BSdat$site_name)
}

# Check for NA values
if (any(is.na(pycno_BSdat$site_name))) {
  pycno_BSdat <- na.omit(pycno_BSdat)
}

# Check unique levels
print(unique(pycno_BSdat$site_name))

# Refit the model
mod_pycno4_glmm_rescaled <- glmer(total ~ comparison_date + kelp.pa + (1 | site_name),
                                  family = poisson,
                                  data = pycno_BSdat)
isSingular(mod_pycno4_glmm_rescaled)
summary(mod_pycno4_glmm_rescaled)




#############GAM
## Q1: Pycnopodia recovery - is there a difference in total between years?
mod_pycno1GAM <-gam(total ~ comparison_date 
                    + s(site_name, bs = "re"),  # Random effect for site_name
                    family = poisson(link = "log"),  # Poisson distribution with log link
                    data = pycno_BSdat,
                    method = "REML") 

summary(mod_pycno1GAM)
AIC(mod_pycno1GAM)
#AIC: 414.3861, p2021 = 0.5408 , p2022=  0.1462   , p2023= 0.0069, 
#high just like gamglss, same relationships 

#when we add kelp presence/absence as an additional predictor
mod_pycno2GAM <- gam(total ~ comparison_date + kelp.pa 
                     + s(site_name, bs = "re"),
                     family = poisson(link = "log"),
                     method = "REML",
                     data = pycno_BSdat)
summary(mod_pycno2GAM)
AIC(mod_pycno2GAM)
#AIC: 409.1768, 
#same relationships as with GAMLSS

#when we add size_class as an additional predictor
mod_pycno3GAM <-gam(total ~ comparison_date + size_class + kelp.pa
                    + s(site_name, bs = "re"),  # Random effect for site_name
                    family = poisson(link = "log"),  # Poisson distribution with log link
                    data = pycno_BSdat,
                    method = "REML")
summary(mod_pycno3GAM)
AIC(mod_pycno3GAM)
#AIC: 332.4017
#same relationships as with GAMLSS
bestGAMmod <- mod_pycno3GAM


######GAMLSS
mod_pycno1 <- gamlss(total ~ comparison_date 
                     + re(random = ~1 | site_name),
                     family = PO(),
                     method = RS(),
                     data = pycno_BSdat, 
                     control = gamlss.control(n.cyc = 200))
summary(mod_pycno1)
#AIC: 414.8407 , p2021 = 0.47016 , p2022=  0.11159, p2023= 0.00339

#when we add kelp presence/absence as an additional predictor
mod_pycno2 <- gamlss(total ~ comparison_date + kelp.pa 
                     + re(random = ~1 | site_name),
                     family = PO(),
                     method = RS(),
                     data = pycno_BSdat, 
                     control = gamlss.control(n.cyc = 200))
summary(mod_pycno2)
#AIC: 410.3823, same relationship across years, kelp.pa = 2.07e-09

#when we add size_class as an additional predictor
mod_pycno3 <- gamlss(total ~ comparison_date + size_class + kelp.pa 
                     + re(random = ~1 | site_name),
                     family = PO(),
                     method = RS(),
                     data = pycno_BSdat, 
                     control = gamlss.control(n.cyc = 200))
summary(mod_pycno3)
#AIC: 336.3757, BEST model, p2021= 0.35151, p2022= 0.80506, p2023 = 0.01039
#all (understory, forest, low) kelp.pa significant 

#No include transect number
mod_pycno4 <- gamlss(total ~ comparison_date + size_class + kelp.pa + transect
                     + re(random = ~1 | site_name),
                     family = PO(),
                     method = RS(),
                     data = pycno_BSdat, 
                     control = gamlss.control(n.cyc = 200))
summary(mod_pycno4)
#transect doesn't have a significant effect, AIC= 338.403


bestmodel <- mod_pycno3
plot(bestmodel)

#We see that significant change in total occurred in 2023 but that there is 
#no significant change in total from early post-SSWD levels (2021, 2022) and 2024









##Q2: Pycnopodia recovery - is there a difference in size between years?
pycno_BSdat$survey_year <- as.factor(pycno_BSdat$survey_year)
pycno_BSdat$size_class <- as.numeric(pycno_BSdat$size_class)
sum(is.infinite(pycno_BSdat$size_class))
sum(is.nan(pycno_BSdat$size_class))
sum(is.infinite(pycno_BSdat$survey_year))
sum(is.nan(pycno_BSdat$survey_year))


aov_result <- aov(size_class ~ survey_year, data=pycno_BSdat)
summary(aov_result)
#Normal?
residuals_aov <- residuals(aov_result)
qqnorm(residuals_aov)
qqline(residuals_aov, col = "red")
shapiro.test(residuals_aov)
# p=1.35e-06; significant deviation from normality

#Do a non-parametric test
KT <- kruskal.test(size_class ~ comparison_date, data = pycno_BSdat)
print(KT)
# distributions of size_class are significantly different across years
# p=0.03311

dunn_results <- dunnTest(size_class ~ comparison_date, data = pycno_BSdat, method = "bonferroni")
print(dunn_results)
#The only significance exists between 2021 and 2023 (p=0.022) where 2021 
#contains organisms from smaller size classes



# Summary data calculation
summary_data <- pycno_BSdat %>%
  group_by(survey_year, size_class) %>%
  summarise(total_count = sum(total, na.rm = TRUE), .groups = 'drop')

# Create a summary for x-axis annotations
total_count_summary <- summary_data %>%
  group_by(survey_year) %>%
  summarise(total_count = sum(total_count), .groups = 'drop')

# Ensure size_class is treated as a factor with specified levels
summary_data <- summary_data %>%
  mutate(size_class = factor(size_class, levels = c(0, 2.5, 5, 7.5, 10, 12.5, 15, 20, 30)))  # Specified order

# Create the plot
size_plot1 <- ggplot(summary_data, aes(x = factor(survey_year), y = size_class)) +
  geom_tile(aes(width = total_count / max(total_count)), color = "black", alpha = 0.7, fill = "lightblue") +
  geom_text(data = total_count_summary, aes(x = factor(survey_year), y = -Inf, label = paste("n=", total_count)), vjust = -0.4, color = "darkgrey") +
  scale_size_continuous(name = "Total Count") +
  theme_classic(base_size = 14) +
  labs(
    x = "Survey Year",
    y = "Size Class (cm)",
    fill = "Survey Year",
    color = "Survey Year"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(size_plot1)















#Plot to show size distribution
summary_data <- pycno_BSdat %>%
  group_by(survey_year, size_class) %>%
  summarise(total_count = sum(total, na.rm = TRUE), .groups = 'drop')

# Create a summary for x-axis annotations
total_count_summary <- summary_data %>%
  group_by(survey_year) %>%
  summarise(total_count = sum(total_count), .groups = 'drop')

View(summary_data)
# Create the plot
size_plot1 <- ggplot(summary_data, aes(x = factor(survey_year), y = size_class)) +
  geom_tile(aes(width = total_count / max(total_count)), color = "black", alpha = 0.7, fill = "lightblue") +
  geom_text(data = total_count_summary, aes(x = factor(survey_year), y = -Inf, label = paste("n=", total_count)), vjust = -0.4, color = "darkgrey") +
  scale_size_continuous(name = "Total Count") +
  theme_classic(base_size = 14) +
  labs(
    x = "Survey Year",
    y = "Size Class (cm)",
    fill = "Survey Year",
    color = "Survey Year"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(size_plot1)










# Summary data calculation, omitting NAs
summary_data <- pycno_BSdat %>%
  filter(!is.na(size_class) & !is.na(total)) %>%  # Omit rows with NAs
  group_by(survey_year, size_class) %>%
  summarise(total_count = sum(total, na.rm = TRUE), .groups = 'drop')

# Create a summary for x-axis annotations
total_count_summary <- summary_data %>%
  group_by(survey_year) %>%
  summarise(total_count = sum(total_count), .groups = 'drop')

# Create the plot
size_plot1 <- ggplot(summary_data, aes(x = factor(survey_year), y = size_class)) +
  geom_tile(aes(width = total_count / max(total_count)), color = "black", alpha = 0.7, fill = "lightblue") +
  geom_text(data = total_count_summary, aes(x = factor(survey_year), y = -Inf, label = paste("n=", total_count)), vjust = -0.4, color = "darkgrey") +
  scale_size_continuous(name = "Total Count") +
  theme_classic(base_size = 14) +
  labs(
    x = "Survey Year",
    y = "Size Class (cm)",
    fill = "Survey Year",
    color = "Survey Year"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(size_plot1)



# Summary data calculation
summary_data <- pycno_BSdat %>%
  group_by(survey_year, size_class) %>%
  summarise(total_count = sum(total, na.rm = TRUE), .groups = 'drop')

# Create a summary for x-axis annotations
total_count_summary <- summary_data %>%
  group_by(survey_year) %>%
  summarise(total_count = sum(total_count), .groups = 'drop')

# Ensure size_class is treated as a factor with levels from 0 to 30
summary_data <- summary_data %>%
  mutate(size_class = factor(size_class, levels = 0:30))  # Adjust levels as needed

# Create the plot
size_plot1 <- ggplot(summary_data, aes(x = factor(survey_year), y = size_class)) +
  geom_tile(aes(width = total_count / max(total_count)), color = "black", alpha = 0.7, fill = "lightblue") +
  geom_text(data = total_count_summary, aes(x = factor(survey_year), y = -Inf, label = paste("n=", total_count)), vjust = -0.4, color = "darkgrey") +
  scale_size_continuous(name = "Total Count") +
  theme_classic(base_size = 14) +
  labs(
    x = "Survey Year",
    y = "Size Class (cm)",
    fill = "Survey Year",
    color = "Survey Year"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(size_plot1)













#Violin plot as initial plot type
summary_data <- pycno_BSdat %>%
  group_by(survey_year, size_class) %>%
  summarise(total_count = sum(total, na.rm = TRUE), .groups = 'drop')

# Create the plot
size_plot2 <- ggplot(summary_data, aes(x = factor(survey_year), y = size_class)) +
  geom_violin(color = "black", alpha = 0.7) +  # Remove fill aesthetic
  scale_size_continuous(name = "Total Count") +
  theme_classic(base_size = 14) +
  labs(
    x = "Survey Year",
    y = "Size Class (cm)"
  ) +
  theme(
    legend.position = "none"  # Remove the legend
  ) 

print(size_plot2)

size_plots <- (size_plot1 + size_plot2) + plot_layout(ncol = 2)
size_plots





#Determine influence of kelp presence on pycnopodia abundance
#merge low kelp.pa with absent
pycno_BSdat <- pycno_BSdat %>%
  mutate(kelp.pa = recode(kelp.pa, "low" = "absent"))

#plot abundance data as column plot, include kelp.pa as a level
pycno_BSdat$density <- pycno_BSdat$total / (pycno_BSdat$transect * 100)

density_plot <- ggplot(pycno_BSdat) +
  aes(fill = kelp.pa, y = density, x = survey_year) +
  theme_classic() +
  geom_bar(position = "stack", stat = "identity") + 
  ylab("Pycnopodia Density (count/100m2)") +
  xlab("Survey Year") + 
  geom_text(aes(x = survey_year, y = 0,  # Position the text below the bars
                label = label_text, group = kelp.pa),
            position = position_dodge(width = 0.9), vjust = 1.0, color = "darkgrey", size = 3.0) +  # Adjust vjust for spacing below bars
  scale_fill_manual(
    values = c("absent" = "#FF6666", "forest" = "#58BC82", "understory" = "#67DD30"),  # Customize these colors
  ) +
  guides(fill = guide_legend(title = "Kelp")) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Pycnopodia Density by Survey Year")

density_plot

# Abundance across kelp & non-kelp sites
# Step 1: Aggregate the data by site first, summing the 'total' at each site

site_level_data <- pycno_BSdat %>%
  group_by(survey_year, kelp.pa, site_name, transect) %>%
  summarize(site_total = sum(density, na.rm = TRUE), .groups = 'drop')  # Ungroup after summarizing

# Calculate the mean density as the mean of the summed site totals
aggregated_data <- site_level_data %>%
  group_by(survey_year, kelp.pa, transect) %>%
  summarize(mean_density = mean(site_total, na.rm = TRUE), .groups = 'drop')  # Ungroup after summarizing

# Apply the division based on kelp.pa
aggregated_data <- aggregated_data %>%
  mutate(mean_density = case_when(
    kelp.pa == "absent" ~ mean_density,
    kelp.pa == "forest" ~ mean_density,
    kelp.pa == "understory" ~ mean_density
  ))

# Step 3: Calculate the number of unique sites for each kelp.pa category
num_sites_data <- site_level_data %>%
  group_by(survey_year, kelp.pa) %>%
  summarize(num_sites = n_distinct(site_name), .groups = 'drop')  # Number of distinct sites
# Step 4: Join the num_sites_data with aggregated_data
aggregated_data <- aggregated_data %>%
  left_join(num_sites_data, by = c("survey_year", "kelp.pa"))

# Create a new label column with "n=" prefix
aggregated_data_filtered <- aggregated_data %>%
  mutate(label_text = paste("n=", num_sites, sep = ""))

# Ensure kelp.pa is a factor with the desired order
aggregated_data_filtered <- aggregated_data_filtered %>%
  mutate(kelp.pa = factor(kelp.pa, levels = c("absent", "understory", "forest")))

# Ensure kelp.pa is a factor with the desired order
aggregated_data_filtered <- aggregated_data_filtered %>%
  mutate(kelp.pa = factor(kelp.pa, levels = c("absent", "understory", "forest")),
         survey_year = factor(survey_year, levels = unique(survey_year)))  # Optional: to ensure survey_year is ordered correctly

# Create the bar plot
densitykpa_plot <- ggplot(data = aggregated_data_filtered) +
  geom_bar(aes(x = survey_year, y = mean_density, fill = kelp.pa), 
           stat = "identity", position = "dodge") +
  geom_text(aes(x = survey_year, y = 0,  # Position the text below the bars
                label = label_text, group = kelp.pa),
            position = position_dodge(width = 0.9), vjust = 1.0, color = "darkgrey", size = 3.0) +  # Adjust vjust for spacing below bars
  labs(x = "Survey Year", y = "Mean Density (count/100m2) per site type", fill = "Kelp Presence") +
  theme_classic() +
  scale_fill_manual(
    values = c("absent" = "#FF6666", "understory" = "#67DD30", "forest" = "#58BC82"),  # Customize these colors
    name = "Macroalgae"  # Set the legend title
  ) +
  theme(
    legend.position = "top",  # Move legend to the top
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Average Density of Pycnopodia at Kelp and Non-Kelp Sites")

print(densitykpa_plot)

View(aggregated_data_filtered)


#Plot the average occupancy of Pycnopodia at Kelp and Non-Kelp Sites
# Step 1: Calculate occupancy at each site
site_level_data <- pycno_BSdat %>%
  group_by(survey_year, kelp.pa, site_name) %>%
  summarize(occupancy = ifelse(sum(total, na.rm = TRUE) > 0, 1, 0), .groups = 'drop')  # 1 for presence, 0 for absence
# Step 2: Calculate the mean occupancy per site
aggregated_data <- site_level_data %>%
  group_by(survey_year, kelp.pa) %>%
  summarize(mean_occupancy = mean(occupancy, na.rm = TRUE), .groups = 'drop')  # Mean of occupancy (proportion of sites with presence)
# Step 3: Calculate the number of unique sites for each kelp.pa category
num_sites_data <- site_level_data %>%
  group_by(survey_year, kelp.pa) %>%
  summarize(num_sites = n_distinct(site_name), .groups = 'drop')  # Number of distinct sites
# Step 4: Join the num_sites_data with aggregated_data
aggregated_data <- aggregated_data %>%
  left_join(num_sites_data, by = c("survey_year", "kelp.pa"))
# Step 5: Filter out rows with missing or zero mean_occupancy
aggregated_data_filtered <- aggregated_data %>%
  filter(!is.na(mean_occupancy))

# Create a new label column with "n=" prefix
aggregated_data_filtered <- aggregated_data_filtered %>%
  mutate(label_text = paste("n=", num_sites, sep = ""))

aggregated_data_filtered <- aggregated_data_filtered %>%
  mutate(kelp.pa = factor(kelp.pa, levels = c("absent", "understory", "forest")),
         survey_year = factor(survey_year, levels = unique(survey_year)))  # Optional: to ensure survey_year is ordered correctly

# Create the plot
occupancy_plot <- ggplot(data = aggregated_data_filtered) +
  geom_bar(aes(x = survey_year, y = mean_occupancy * 100, fill = kelp.pa), 
           stat = "identity", position = "dodge") +
  geom_text(aes(x = survey_year, y = 0,  # Position the text below the bars
                label = label_text, group = kelp.pa),
            position = position_dodge(width = 0.9), vjust = 1.0, color = "darkgrey", size = 3.0) +  # Adjust vjust for spacing below bars
  labs(x = "Survey Year", y = "Site Pycnopodia Observed (%)") +
  theme_classic() +
  scale_fill_manual(
    values = c("absent" = "#FF6666", "understory" = "#67DD30", "forest" = "#58BC82"),  # Customize these colors
    name = "Macroalgae"  # Set the legend title
  ) +
  theme(
    legend.position = "top",  # Move legend to the top
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Average Occupancy of Pycnopodia at Kelp and Non-Kelp Sites")

# Display the plot
print(occupancy_plot)

kelp_plots <- (densitykpa_plot + occupancy_plot) + plot_layout(ncol = 2)
kelp_plots






#Occupancy vs mean total
#Every point represents all sites averaged over 2021-2024 for each kelp presence category

# Step 1: Calculate summed total per site and then the mean of these sums by site and kelp presence/absence
summed_total_by_site <- pycno_BSdat%>%
  group_by(site_name) %>%
  summarise(
    summed_total = sum(total, na.rm = TRUE),
    .groups = "drop"
  )

# Now calculate the mean of these summed totals by kelp presence/absence
aggregated_data <- pycno_BSdat %>%
  left_join(summed_total_by_site, by = "site_name") %>%
  group_by(kelp.pa) %>%
  summarise(
    occupancy = mean(total > 0),  # Proportion of Pycnopodia presence per site
    mean_summed_total = mean(summed_total),  # Mean of summed totals per site group
    .groups = "drop"
  )

# Step 2: Create the plot
OvsA_plotKPA <- ggplot(aggregated_data, aes(y = occupancy*100, x = mean_summed_total, color = kelp.pa)) +
  geom_point(size = 6, alpha = 1.2) +
  scale_color_manual(values = c("absent" = "#FF6666", "understory" = "#67DD30", "forest" = "#58BC82")) +  # Custom colors for kelp
  labs(
    y = "Percentage of Sites Occupied (%)",
    x = "Mean Total Abundance per site",
    color = "Kelp"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ylim(0, 100.00) +
  ggtitle("Occupancy vs. Mean total at Kelp and Non-Kelp Sites")

# Print the plot
OvsA_plotKPA




#Occupancy vs mean total
#Every point represents the summed total abundance/occupancy across sites of 
#each kelp presence category for each year

# Step 1: Calculate summed total per site and then the mean of these sums by site and kelp presence/absence
summed_total_by_site <- pycno_BSdat%>%
  group_by(site_name) %>%
  summarise(
    summed_total = sum(total, na.rm = TRUE),
    .groups = "drop"
  )

# Now calculate the mean of these summed totals by kelp presence/absence
aggregated_data <- pycno_BSdat %>%
  left_join(summed_total_by_site, by = "site_name") %>%
  group_by(kelp.pa, survey_year) %>%
  summarise(
    occupancy = mean(total > 0),  # Proportion of Pycnopodia presence per site
    mean_summed_total = mean(summed_total),  # Mean of summed totals per site group
    .groups = "drop"
  )

# Step 2: Create the plot
OvsA_plotYEAR <- ggplot(aggregated_data, aes(y = occupancy * 100, x = mean_summed_total, color = kelp.pa, shape = survey_year)) +
  geom_point(size = 6, alpha = 0.7) +
  scale_color_manual(values = c("absent" = "#FF6666", "understory" = "#67DD30", "forest" = "#58BC82")) +  # Custom colors for kelp
  scale_shape_manual(values = c(16, 17, 15, 18)) +  # Circle, Triangle, Square, Diamond
  labs(
    y = "Percentage of Sites Occupied (%)",
    x = "Mean Total Abundance per site",
    color = "Kelp",
    shape = "Survey Year"  # Add shape legend title
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ylim(0, 100.00) +
  ggtitle("Occupancy vs. Mean Total at Kelp and Non-Kelp Sites")

print(OvsA_plotYEAR)



#size indicaing n=sites
# Step 1: Calculate summed total per site and then the mean of these sums by site and kelp presence/absence
summed_total_by_site <- pycno_BSdat %>%
  group_by(site_name) %>%
  summarise(
    summed_total = sum(total, na.rm = TRUE),
    .groups = "drop"
  )

# Now calculate the mean of these summed totals by kelp presence/absence and count unique site_name
aggregated_data <- pycno_BSdat %>%
  left_join(summed_total_by_site, by = "site_name") %>%
  group_by(kelp.pa, survey_year) %>%
  summarise(
    occupancy = mean(total > 0),  # Proportion of Pycnopodia presence per site
    mean_summed_total = mean(summed_total),  # Mean of summed totals per site group
    unique_site_count = n_distinct(site_name),  # Count of unique site_name
    .groups = "drop"
  )

# Step 3: Create a size category based on unique site counts
aggregated_data <- aggregated_data %>%
  mutate(size_category = case_when(
    unique_site_count < 6 ~ "1-5",   
    unique_site_count < 19 ~ "6-18",     
    TRUE ~ "18+"
  ))


# Create the plot with continuous size scaling
OvsA_plot <- ggplot(aggregated_data, aes(y = occupancy * 100, x = mean_summed_total, 
                                         color = kelp.pa, shape = survey_year, size = size_category)) +
  geom_point(alpha = 0.7) +  
  scale_color_manual(values = c("absent" = "#FF6666", "understory" = "#67DD30", "forest" = "#58BC82")) +  
  scale_shape_manual(values = c(16, 17, 15, 18)) +  
  scale_size_manual(values = c(`1-5` = 4, `6-18` = 6, `18+` = 10)) +  # Correctly define sizes for each category
  labs(
    y = "Percentage of Sites Occupied (%)",
    x = "Mean Total Abundance per site",
    color = "Kelp",
    shape = "Survey Year",  
    size = "Number of Sites"  
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  ylim(0, 100.00) +
  ggtitle("Occupancy vs. Mean Total at Kelp and Non-Kelp Sites") +
  guides(color = guide_legend(override.aes = list(size = 5)),  # Adjust the size of points in the color legend
         shape = guide_legend(override.aes = list(size = 5)))   # Adjust the size of shapes in the shape legend
         

print(OvsA_plot)

#make a panel
# Combine the plots with patchwork
top_row <- (densitykpa_plot + occupancy_plot) + plot_layout(ncol = 2)

topL <-top_row + plot_annotation(
  tag_levels = 'A',  # Set to 'A' to use labels A, B, C, etc.
  tag_suffix = ')'  # Adds a colon after the label (e.g., A:, B:, C:)
) 


OvsA_plot



# Size distribution VI
model_pycno_VISIZE <- pycno_dat %>%
  filter(survey_year == "2023")

# Tile plot to show size distribution, excluding "no size"
summary_dataVIS <- model_pycno_VISIZE %>%
  filter(size_class != "no size") %>%
  group_by(survey_region, size_class) %>%
  summarise(total_count = sum(total, na.rm = TRUE), .groups = 'drop')

# Rename survey regions
summary_dataVIS <- summary_dataVIS %>%
  mutate(survey_region = recode(survey_region,
                                "BS" = "Barkley Sound",
                                "CR" = "Georgia Strait",
                                "RR" = "Race Rocks",
                                "Vic" = "Victoria"))

# Set the order of size_class on the y-axis
size_order <- c("0", "2.5", "5", "7.5", "10", "12.5", "15", "20", "30", "35", "40", "50")
summary_dataVIS$size_class <- factor(summary_dataVIS$size_class, levels = size_order)

# Create a summary for x-axis annotations
total_count_summaryVIS <- summary_dataVIS %>%
  group_by(survey_region) %>%
  summarise(total_count = sum(total_count), .groups = 'drop')

# Create the plot
size_plot3 <- ggplot(summary_dataVIS, aes(x = survey_region, y = size_class)) +
  geom_tile(aes(width = total_count / max(total_count)), color = "black", alpha = 0.7, fill = "lightblue") +
  geom_text(data = total_count_summaryVIS, aes(x = survey_region, y = -Inf, label = paste("n=", total_count)), vjust = -1, color = "darkgrey") +
  scale_size_continuous(name = "Total Count") +
  theme_classic(base_size = 14) +
  labs(
    x = "Survey Region",
    y = "Size Class (cm)",
    fill = "Survey Region"
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)) 
  

print(size_plot3)



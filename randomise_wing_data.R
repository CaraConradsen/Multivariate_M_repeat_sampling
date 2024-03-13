library(data.table); library(dplyr)

wing_dat <- fread("mr_wings_6traits.csv")

summary(wing_dat)

unique(wing_dat$Trait)

for (i in unique(wing_dat$Trait)) {
 print(mean(wing_dat[Trait==i, Score]))
}

#Is generation six (Jack's Generation) different?
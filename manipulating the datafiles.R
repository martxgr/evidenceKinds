library(tidyverse)

#study 1
#raw qualtrics csv
d <- read.csv("lb.s1.csv")
#remove metadata rows
#only good pids
d <- d[d$prolific %in% d$pid, ]
d <- d %>%
  dplyr::select(-pid)
#replacing participant id row with index
d$prolific <- seq_len(nrow(d))
#write csv
write.csv(d, "Study1.csv", row.names = FALSE)

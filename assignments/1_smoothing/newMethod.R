library(tidyverse)
library(splines)
library(latex2exp)
library(microbenchmark)

knitr::opts_chunk$set(cache = TRUE, dev.args = list(bg = 'transparent'), 
                      fig.align = "center", fig.pos="h", cache.lazy = TRUE,
                      out.width = "70%")
theme_set(theme_bw(base_size=12))

Nuuk <- read_table2("nuuk.dat.txt", 
                    col_names = c("Year", 1:12), 
                    na = "-999", 
                    skip = 1) %>% 
  gather(key = "Month", value = "Temperature", -Year, convert = TRUE) %>% 
  mutate(Temperature = Temperature / 10) %>% 
  filter(Year > 1866)
Nuuk_year <- group_by(Nuuk, Year) %>% 
  summarise(
    Temperature = mean(Temperature),
    Median = median(Temperature),
    High = max(Temperature), 
    Low = min(Temperature)
  )
n <- nrow(Nuuk_year)



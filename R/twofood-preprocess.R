setwd("~/Dropbox/research/lsjm-art/lsjm-code")

library(dplyr)
library(magrittr)

source("R/art-functions.R")

## fix_num and roi vary per item and person -> it might mean something, but I drop them (used distinct without these columns)
load("data/SmithKrajbich2018.RData")

df = data %>% select(SubjectNumber, LeftRight, RT, ValueLeft, ValueRight)
df = distinct(df)
names(df) = c("person", "item", "resp", "RT" )
name_item = unique(df$item)
name_person = unique(df$person)
pick_item = df$item %in% name_item

df = df[pick_item, ]

pick_person = df$person %in% unique(df$person)
df = df[pick_person,]

di = df[,-4]
dt = df[,-3]

nitem = length(unique(df$item))
nperson = length(unique(df$person))

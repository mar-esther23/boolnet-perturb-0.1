# [*] simplify getFixedAttractors
# [*] generalize attractor2dataframe
# [*] rewrite fixedAttractorsListsToDataframe

require(BoolNet)
source("BoolNet-extensions.R")

net <- loadNetwork("minTh17iTreg.txt")
attr <- getAttractors(net)

labels.rules <- data.frame( labels = c('Th0', 'Th17', 'Treg', 'IL10+', 'TGFB+', 'RORGT+'),   rules  = c('!(RORGT | FOXP3 | TGFB | IL10)',  'RORGT & STAT3',  'FOXP3 & TGFB',  'IL10',  'TGFB & ! (RORGT | FOXP3)',  'RORGT & ! STAT3' ), stringsAsFactors = FALSE)
labels <- labelAttractors(attr, net$genes, labels.rules$labels, labels.rules$rules)
attr.labels <- attr

for (i in 1:length(attr.labels$attractors))  attr.labels$attractors[[i]]$label <- labels[[i]]










mutants <- getFixedAttractors(net)
mutants

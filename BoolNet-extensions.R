

dec2binState <- function(x, genes){ 
    # Takes an integer and a list of node names
    # Returns a vector of 0s and 1s where the name of elemento corresponds to a named node
    state <- as.integer( intToBits(x)[1:length(genes)] )  
    names(state) <- genes
    state
}

######################
####   LABELING   ####
######################
labelState <- function(state, node.names, labels, rules) {
    # label a single binary state
    # returns a label string
    names(state) <- node.names
    label = c()
    for (j in 1:length(rules)) { #evaluate rules
        # create string with function
        f <- paste('function(', paste(node.names, collapse=','), ') { if (',
                   rules[j], ') \'', labels[j], '\' }' , sep='')
        f <- eval(parse(text=f)) # evaluate string to create function
        label <- append(label, do.call(f, as.list(state))) # apply function
    }
    # format label
    if (is.null(label)) { label <-c('NoLabel')}
    label <- paste(label, collapse='')
}

labelAttractors <- function(attr, node.names, labels, rules) {
    # takes an attractors object created by BoolNet
    # returns a list of the labels for each attractor in order.
    # If an attractor has multiple states it will return a label for each state.
    res <- list()
    for (i in 1:length(attr$attractors)) {
        label <- sapply(attr$attractors[[i]]$involvedStates, function(state) {
            state <- dec2binState(state, node.names) #state to binary
            l <- labelState(state, node.names, labels, rules) #label
        })
        res <- append(res, list(label))
    }
    res
}
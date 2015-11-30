

dec2binState <- function(x, genes){ 
    # Takes an integer and a list of node names
    # Returns a vector of 0s and 1s where the name of elemento corresponds to a named node
    state <- as.integer( intToBits(x)[1:length(genes)] )  
    names(state) <- genes
    state
}

attractor2dataframe <- function(attr) {
    # Convert an BoolNet attractor object to a data frame.
    # states will be converted to strings and collapsed using '/'
    data.frame(
        involvedStates = sapply(attr$attractors, function(a) {
            paste(as.character(a$involvedStates), collapse='/')
        })  , 
        basinSize = sapply(attr$attractors, function(a) a$basinSize )
    )}



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



#######################
####   FIX NODES   ####
#######################
getStatesAndBasin <- function (attr) {
    # Receives a set of attractors, 
    # returns a list with int_state and basin
    # if the basin value is NA returns TRUE
    states <- sapply(attr$attractors, function(attractor) attractor$involvedStates )
    basin <- sapply(attr$attractors, function(attractor) {
        if (is.na(attractor$basinSize)) TRUE  #control for asynchronous
        else attractor$basinSize  #basin size in synchronous
    })
    states <- list(states=states, basin=basin)
    #names(states) = states
    states
}

getFixedAttractors <- function(net, genes, value, label, type="synchronous", returnDataFrame=TRUE) {
    #     Simulates knocked-out or over-expression for all  genes
    #     by fixing the values of genes to 0 or 1, 
    #     Returns all the possible attractors with basin if posible
    
    # generate genes to evaluate
    if (missing(genes) | missing(value)) { #Default, evaluate all single KOver
        genes <- c(NA, net$genes, net$genes)
        value <- c(NA, rep(0, length(net$genes)), rep(1, length(net$genes)))
    }
    if (missing(label)) {
        label <- paste(genes, value, sep='_')
        label[match('NA_NA', label)] <- 'WT'
    }
    
    # Calculate mutants
    mutants <- list()
    for (i in 1:length(genes)) {
        #print(paste(i, label[i], genes[i], value[i]  ))
        if (!is.na(genes[i])) { net <- fixGenes(net, unlist(genes[i]), unlist(value[i])) }
        attr <- getAttractors(net, type=type)
        states <- getStatesAndBasin(attr)
        mutants[[i]] <- states
        if (! is.na(unlist(genes[i]))) net <- fixGenes(net, unlist(genes[i]), -1)
    }
    # convert to dataframe
    names(mutants) <- label
    if (returnDataFrame==TRUE) mutants <- fixedAttractorsListsToDataframe(mutants)
    mutants
}

fixedAttractorsListsToDataframe <- function(attrList) {
    # Receives a list of lists
    # Each list corresponds to a network (with != fixed genes)
    # Each list contains lists $states and $basin
    # $states can cointain cyclic attractors
    # Return a dataframe
    
    # Convert cycles to str
    for (i in 1:length(attrList)) {
        for (j in 1:length(attrList[[i]][[1]])) {
            s <- attrList[[i]][[1]][[j]]
            if (length(s) > 1) { s <- paste( s, collapse='-') }
            attrList[[i]][[1]][[j]] <- s }
        attrList[[i]][[1]] <- unlist(attrList[[i]][[1]])
    }
    
    # Determine dataframe dimensions
    states <- sapply(attrList, function(attr) attr$states)
    states <- sort(unique(unlist(states)))
    df <- data.frame(attr=states, row.names=states, stringsAsFactors = FALSE)
    
    #merge
    for (name in names(attrList)) {
        attr <- data.frame(attrList[[name]], row.names=attrList[[name]]$states, stringsAsFactors = FALSE)
        names(attr) <- c('attr', name)
        df <- merge(df, attr, all=TRUE)
    }
    df$attr <- as.character(df$attr)
    row.names <- as.character(df$attr)
    df
}
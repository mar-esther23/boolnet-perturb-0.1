

dec2binState <- function(x, genes){ 
    # Takes an integer and a list of node names
    # Returns a vector of 0s and 1s where the name of elemento corresponds to a named node
    state <- as.integer( intToBits(x)[1:length(genes)] )  
    names(state) <- genes
    state
}



attractor2dataframe <- function(attr) {
    # Convert an BoolNet attractor object to a data frame of attr$attractor properties.
    # attr$attractor properties with multiple elements will be transformed to strings and joined with "/"
    
    attr <- attr$attractors
    # create properties list, if labeled we will have more
    attr.properties <- vector("list", length(attr[[(1)]]))
    names(attr.properties) <- names(attr[[(1)]])
    attr.properties
    
    for (n in names(attr.properties) ) { #create list for each property
        attr.properties[[n]] <- sapply(attr, function(a) a[[n]]) 
        #verify number of elements inside list
        ncol <- max(sapply(attr.properties[[n]], length))
        if ( ncol > 1) { #collapse
            attr.properties[[n]] <- sapply(attr.properties[[n]], function(a) {
                paste(as.character(a), collapse='/')
            })}    
    }
    data.frame(attr.properties)
}



attractorsListsToDataframe <- function(attr.list) {
    # Receives a list of BoolNet attractors and return a dataframe
    # Each column is named attrName.propertyName
    
    # Transform from attr to dataframe
    attr.list <- lapply(attr.list, attractor2dataframe)
    # set involvedStates as rowname, delete
    # and rename df columns to attrName.propertyName
    for (n in names(attr.list)) {
        rownames(attr.list[[n]]) <- attr.list[[n]]$involvedStates #set involvedStates as rowname
        attr.list[[n]]$involvedStates <- NULL #delete
        names(attr.list[[n]]) <- paste(n, names(attr.list[[n]]), sep='.') # rename df columns
    } 
    
    #merge and reduce by rownames
    attr.df <- Reduce(function(x, y){
        df <- merge(x, y, by= "row.names", all=TRUE)
        rownames(df) <- df$Row.names
        df$Row.names <- NULL
        return(df)
    }, attr.list)
    attr.df
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



#######################
####   FIX NODES   ####
#######################
getFixedAttractors <- function(net, genes, value, label, type="synchronous", returnDataFrame=TRUE) {
    # Takes a net and fixes the genes with value, returns attractors
    # net:      network
    # genes:    list of genes to fix
    # values:   list of values to fix genes
    # label:    names of fixed networks
    # type:     update type, if async the basins are T of F
    # returnDataFrame: datatype to return, if true dataframe, if false list of attractors
    
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
        mutants[[i]] <- getAttractors(net, type=type)
        if (! is.na(unlist(genes[i]))) net <- fixGenes(net, unlist(genes[i]), -1)
    }
    names(mutants) <- label
    # convert attractors to dataframe
    if (returnDataFrame==TRUE) mutants <- attractorsListsToDataframe(mutants)
    mutants
}




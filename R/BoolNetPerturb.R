#' Convert integer to binary state vector with node names.
#' '
#' @param x input integer representing the state
#' @param node.names network node names
#' @return Numeric binary vector of the same length as the number of nodes. Each 
#'     position corresponds to a node of the network. The values of each element 
#'     are 0 or 1. The name of each element corresponds to the name of the node in
#'     that position in the network.
#' @seealso \code{\link{intToBits}} which this function wraps
#' @export
#' @examples
#' int2binState(0, net$genes)
int2binState <- function(x, node.names){ 
    state <- as.integer( intToBits(x)[1:length(node.names)] )  
    names(state) <- node.names
    state
}



#' Determine if a node is an input of the network.
#'
#' @param gene node to asses
#' @param net BoolNet network 
#' @return Boolean is the node an input in the network.
#' @export
#' @examples
#' isGeneInput("cycd", cellcycle)
isGeneInput <- function(gene, net) {
    if ( 
        all((which( net$genes == gene )) == net$interactions[[gene]]$input)
        && 
            all(c(0,1) == net$interactions[[gene]]$func)
    ) { return(TRUE) } else return(FALSE)
}



#' Modify state, set new value to target nodes.
#'
#' @param state original state
#' @param new.nodes nodes to modify
#' @param new.values new values of modified nodes
#' @return State with new values in modified nodes.
#' @export
#' @examples
#' state <- int2binState(0, cellcycle$genes)
#' new.state <- setStateValues(state, c("cycd","cyce"), c(0,1))
setStateValues <- function(state, new.nodes, new.values) {
    for (i in 1:length(new.nodes)) { state[new.nodes[i]] <- new.values[i] }
    state
}



####################
#### DATAFRAMES ####
####################

#' Convert attractor to data frame.
#' 
#' Convert a BoolNet attractor object to a data frame. Each property of attr$attractors corresponds to a dataframe column. If the property has elements with length > 1 it converts them to a string and joins them with sep.
#'
#' @param attr BoolNet attractor object
#' @param sep string to join elements with length > 1, default "/"
#' @return Dataframe, each column corresponds to a property of the attractor
#' @export
#' @examples
#' > attr <- getAttractors(cellcycle$genes)
#' > attractor2dataframe(attr)
#'               involvedStates basinSize
#' 1                        162       512
#' 2 25/785/849/449/389/141/157       512
attractor2dataframe <- function(attr, sep="/") {
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
                paste(as.character(a), collapse=sep)
            })}    
    }
    data.frame(attr.properties, stringsAsFactors=FALSE)
}



#' Convert a list of attractors to a data frame.
#' 
#' Convert a list of BoolNet attractor objects to a data frame. Each property of eacn attr$attractors corresponds to a dataframe column. Columns are named attrName.propertyName, if the list has no names numbers will be used. If the property has elements with length > 1 it converts them to a string and joins them with sep.
#'
#' @param attr.list list of BoolNet attractor objects
#' @param sep string to join elements with length > 1, default "/"
#' @return Dataframe, each column corresponds to a property of the attractor
#' @export
#' @examples
#' attractorsLists2Dataframe(attr.list)
attractorsLists2Dataframe <- function(attr.list, sep='/') {
    # Receives a list of BoolNet attractors and return a dataframe
    # Each column is named attrName.propertyName
    
    # Transform from attr to dataframe
    attr.list <- lapply(attr.list, attractor2dataframe, sep)
    # Verify names exist
    if ( is.null(names(attr.list)) ) names(attr.list) <- 1:length(attr.list)
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

#' Label a state.
#' 
#' Labels a binary state using a set of labelling rules.
#' The elements of the state correspond the network node names. If a rule is satisfied the corresponding label is appended.
#'
#' @param state binary state to label
#' @param node.names node names of the state, the length must be the same that the state's
#' @param labels label that will be appended if the corresponding rule is TRUE
#' @param rules rules used for labelling the state. If the state satisfies the rule the corresponing label will be appended.  All the node names present in the rules must be in node.names
#' @param sep string to separate the labels when more than one can be applied to the state.
#' @return String corresponding to the label of the state.
#' @export
#' @examples
#' state <- int2binState(0, net$genes)
#' labelState(state, net$genes, labels, rules, sep='')
labelState <- function(state, node.names, labels, rules, sep='') {
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
    label <- paste(label, collapse=sep)
}



#' Label attractor.
#' 
#' Returns a list of labels for a the involved states in a BoolNet attractor object using a set of labelling rules, returns a list of labels.
#'
#' @param attr Boolnet attractor object
#' @param node.names node names of the network
#' @param labels label that will be appended if the corresponding rule is TRUE
#' @param rules rules used for labelling the state. If the state satisfies the rule the corresponing label will be appended.  All the node names present in the rules must be in node.names
#' @param sep string to separate the labels when more than one can be applied to the state.
#' @return List of strings corresponding to the label of the attractor, if an attractor has multiple states it will return a list of strings for that state.
#' @seealso \code{\link{labelState}} which this function wraps
#' @export
#' @examples
#' attr <- getAttractors(net$genes)
#' labelAttractors(attr, net$genes, labels, rules, sep='')
labelAttractors <- function(attr, node.names, labels, rules, sep='') {
    # takes an attractors object created by BoolNet
    # returns a list of the labels for each attractor in order.
    # If an attractor has multiple states it will return a label for each state.
    res <- list()
    for (i in 1:length(attr$attractors)) {
        label <- sapply(attr$attractors[[i]]$involvedStates, function(state) {
            state <- int2binState(state, node.names) #state to binary
            l <- labelState(state, node.names, labels, rules, sep='') #label
        })
        res <- append(res, list(label))
    }
    res
}



#############################
####   PERTURB NETWORK   ####
#############################

#' Perturb network functions.
#' 
#' Simulates fixed function perturbations (knock-out and over-expression) in a network. Takes a network, a list of gene names and corresponding values to perturb, fixes those values in each network and returns the corresponding attractors for all perturbed networks. By default it  returns all the single node knock-out and over-expression of a network.
#'
#' @param net network to perturb
#' @param genes list of gene names to perturb. To perturb multiple nodes at the same time use a vector inside the list.
#' @param value list of values of the perturbed genes. Knock-out is 0, over-expression is 1.
#' @param label name of the perturbation
#' @param type update type, can be "synchronous" (default) or "asynchronous"
#' @param returnDataFrame if TRUE returns a dataframe where the rownames correspond to the states and columns correspond to the attractor basin size of the perturbed network, if FALSE returns a list of BoolNet attractors
#' @return dataframe or list of attractors of the perturbed networks
#' @seealso \code{\link{fixGenes}} 
#' @export
#' @examples
#' # All single gene knock-out and over-expression of a network
#' perturbNetworkFixedNodes(net)
#' 
#' # Cellcycle network CycD over-expression
#' perturbNetworkFixedNodes( cellcycle, list("CycD"), list(1) )
#' 
#' # Cellcycle network CycDknock-out with Rb over-expression
#' perturbNetworkFixedNodes(  cellcycle, list(c("CycD", "Rb")), list(c(0,1)) , list("CycD-RB+") )
#' 
#' # All CycD double gene knock-outs
#' double.KO <- lapply(cellcycle$genes[-1], function(gene) c("CycD", gene))
#' number.double.KO <- length(double.KO)
#' perturbNetworkFixedNodes(  net=cellcycle, genes=double.KO, 
#'                         value=as.list(rep(0, number.double.KO))  ) 
perturbNetworkFixedNodes <- function(net, genes, value, label, type="synchronous", returnDataFrame=TRUE) {
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
#         print(i)
#         print(label[i])
#         print(genes[i])
#         print(value[i] )
        not.WT <- ! (is.null(genes[i]) || is.na(genes[i]))
#         print(not.WT)
        if ( not.WT ) { net <- fixGenes(net, unlist(genes[i]), unlist(value[i])) }
        mutants[[i]] <- getAttractors(net, type=type)
        if ( not.WT ) net <- fixGenes(net, unlist(genes[i]), -1)
    }
    names(mutants) <- label
    # convert attractors to dataframe
    if (returnDataFrame==TRUE) mutants <- attractorsLists2Dataframe(mutants)
    mutants
}



##########################
####   PERTURB PATH   ####
##########################

#' Perturb network state or trajectory.
#' 
#' Modifies the values of a state and returns the resulting attractor or trajectory. If time=NULL the perturbation will be fixed permanently, if time=n the perturbation will last n time steps and then the rules will return to their original values.
#'
#' @param state binary state to perturb
#' @param net BoolNet network to perturb
#' @param genes genes to perturb
#' @param values value of perturbed genes
#' @param time time of perturbation. If time=NULL the perturbation will be fixed permanently, if time=n the perturbation will last n time steps and then the rules will return to their original values.
#' @param returnTrajectory if TRUE returns the trajectory of state perturbation, if FALSE returns the reached attractor. Default is FALSE
#' @return If returnTrajectory=TRUE returns the trajectory of the state perturbation. If FALSE returns the attractor the state reached. Al the updates are synchronous.
#' @seealso \code{\link{...}} 
#' @export
#' @examples
#' # Perturb an state by permanetly fixing the value of nodes, return the final attractor
#' state <- int2binState(162, cellcycle$genes)
#' new.attr <- perturbPathToAttractor(state, cellcycle, 
#'                  c('CycD', 'Rb'), c(1,0))
#' 
#' 
#' # Perturb an state by permanetly fixing the value of nodes, return the trajectory
#' state <- int2binState(162, cellcycle$genes)
#' new.traj <- perturbPathToAttractor(state, cellcycle, 
#'                  c('CycD', 'Rb'), c(1,0), 
#'                  returnTrajectory=TRUE)
#' 
#' # Perturb an state by fixing the value of nodes for one time step, return the final attractor
#' state <- int2binState(162, cellcycle$genes)
#' new.attr <- perturbPathToAttractor(state, cellcycle, 
#'                  c('CycD', 'Rb'), c(1,0), time=1)
#' 
#' 
#' # Perturb an state by fixing the value of nodes one time step, return the trajectory
#' state <- int2binState(162, cellcycle$genes)
#' new.traj <- perturbPathToAttractor(state, cellcycle, 
#'                  c('CycD', 'Rb'), c(1,0), time=1, 
#'                  returnTrajectory=TRUE)
#' 
perturbPathToAttractor <- function(state, net, genes, values, time=NULL, returnTrajectory = FALSE) {    
    initial.state <- state #add conversion later
    names(initial.state) <- net$genes
    if (returnTrajectory) path.perturbed <- list(initial.state) #save initial state in path
    
    if (!is.null(time)) { # if transient perturbation
        # determine original value of inputs
        inputs <- unlist(sapply(net$genes, function(g) {
            if (isGeneInput(g,net)) return(T) else F
        }))
        # We need to save the inputs, because the default input function is input=input
        # So if we don't explicitly return them to the original value 
        # they will stay in the modified value !
        inputs.values <- state[inputs]
        inputs <- names(inputs.values)
        
        # apply perturbation
        net <- fixGenes(net, genes, values) 
        initial.state <- setStateValues(state, genes, values)
        
        for (t in 1:time) { #iterate n times
            if (returnTrajectory) path.perturbed <- append(path.perturbed, list(initial.state))
            initial.state <- stateTransition(net, initial.state)
        }
        # recover original network and inputs
        net <- fixGenes(net, genes, -1) 
        initial.state <- setStateValues(initial.state, inputs, inputs.values) 
    } else { # if fixed perturbation
        # apply perturbation
        net <- fixGenes(net, genes, values) 
        # change the values according to perturbation
        initial.state <- setStateValues(state, genes, values)
        #if (returnTrajectory) path.perturbed <- append(path.perturbed, list(initial.state))
    }
    
    if (!returnTrajectory) { #just return the attractor
        attr <- getAttractors(net, startStates = list(initial.state), returnTable=F)
        return(attr)
    }
    
    # calculate the rest of the trajectory
    path <- getPathToAttractor(net, initial.state)
    
    #if (is.null(time)) return(path)  # if fixed perturbation and trajectory
    
    # join both paths
    path.perturbed <- sapply(path.perturbed, function(f) f) # simplify path.perturbed
    path <- t(path)
    path.res <- t(cbind(path.perturbed, path)) #join lists
    rownames(path.res) <- c(1:(dim(path.res)[1]))
    return(path.res)
}

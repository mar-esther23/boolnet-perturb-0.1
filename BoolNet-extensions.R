

dec2binState <- function(x, genes){ 
    state <- as.integer( intToBits(x)[1:length(genes)] )  
    names(state) <- genes
    state
}

labelState <- function(state, labels, rules) {
    # Evaluate states with rules
    label = c()
    for (j in 1:length(rules)) { #evaluate rules
        # this creates a string with a function using rules-str
        f <- paste('function(', paste(genes, collapse=','), ') { if (', rules[j], ') \'', labels[j], '\' }' , sep='')
        #evaluate string to create function
        f <- eval(parse(text=f))
        #apply function to label
        label <- append(label, do.call(f, s))
    }}
    # format label
    if (is.null(label)) { label <-c('NoLabel') }
    label <- paste(label, collapse='/')
}
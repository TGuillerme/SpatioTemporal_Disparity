##########################
#S3 method for plotting results of various functions of std
##########################
#Plotting the pco load, the global pco or the po through time
#v0.1
##########################
#SYNTAX :
#<data> data from a std function, must be of class "pcoa" for plot.load, "pco.scores" for plot.pco and "pco.slice" for plot.slice
#<legend> whether to add the legend on the plot (default=FALSE)
#<pos.leg> a vector of two numerical values of where to put the legend box. If missing, the legend is plotted by default in the upper left corner. Is ignored if legend=FALSE
#<xlim> & <ylim> two vectors of two values each for setting the limits of the graph axis. If missing, the axis are set to be the maximum/minimum of the data along the axes.
#<col> a list of colors for the plots. If missing, colours are taken from the default palette
#<pars> window parameters for the various plots layout
#<...> any optional argument to be passed to the plot function.
##########################
#----
#guillert(at)tcd.ie 09/10/2014
##########################

plot.std<-function(data, legend=FALSE, pos.leg, xlim, ylim, col, pars=NULL, ...) {

    #SANITIZING
    #data
    class.load<-FALSE
    class.pco<-FALSE
    class.slice<-FALSE
    if(class(data) == 'pcoa') {
        class.load<-TRUE
    } else {
        if(class(data) == 'pco.scores') {
            class.pco<-TRUE
        } else {
            if (class(data) == 'std.slices') {
                class.slice<-TRUE
            } else {
                stop("Data must be from class \"pcoa\", \"pco.scores\" or \"std.slices\"")
            }
        }
    }

    #legend
    check.class(legend, 'logical', ' must be logical.')

    #pos.leg
    if(legend == TRUE) {
        if(missing(pos.leg)) {
            pos.leg='default'
        } else {
            check.class(pos.leg, 'numeric', ' must be a vector of two numerical values.')
            check.length(pos.leg, 2, ' must be a vector of two numerical values.')
        }
    } else {
        pos.leg='none'
    }

    #xlim
    if(missing(xlim)) {
        xlim='default'
    } else {
        check.class(xlim, 'numeric', ' must be a vector of two numerical values.')
        check.length(xlim, 2, ' must be a vector of two numerical values.')    
    }

    #ylim
    if(missing(ylim)) {
        ylim='default'
    } else {
        check.class(ylim, 'numeric', ' must be a vector of two numerical values.')
        check.length(ylim, 2, ' must be a vector of two numerical values.')
    }

    #col
    if(missing(col)) {
        palette("default")
        col<-palette()
    } else {
        check.class(col, 'character', ' col must be a vector of colours.')
    }

    #PLOTTING THE DATA

    if(class.load == TRUE) {
        plot.load(data, legend, pos.leg, ...)
    }

    if(class.pco == TRUE) {
        plot.pco(data, legend, pos.leg, xlim, ylim, col, ...)
    }

    if(class.slice == TRUE) {
        pco.slice(data, legend, pos.leg, xlim, ylim, col, pars, ...)
    }

#End
}
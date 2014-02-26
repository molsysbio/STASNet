## Two functions that simplify data parsing and model creation and fitting

# Global variable to have more outputs
verbose = TRUE;
debug = TRUE;

create_model <- function(model.links="links", data.stimulation="data", basal_activity = "basal.dat")
{
# Creates a parametrized model from an experiment file and the network structure
# It requires the file network_reverse_engineering-X.X/r_binding/fitmodel/R/generate_model.R of the fitmodel package
# model.links should be an adjacency list file representing the network
# Experiment file data.simulation syntax should be as follows, with one line per replicate
#          stimulator                |          inhibitor                |                 	    type                       | [one column per measured nodes]
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
# stimulator name or solvant if none | inhibitor name or solvant if none | c for control, b for blank and t for experiment |    measure for the condition

	### READ DATA
    print("Reading data")
	# Creation of the structure object
	links = read.delim(model.links, header=FALSE)
	model.structure=getModelStructure(links)

	# Read the experiment design and extract the values
	data.file = read.delim(data.stimulation)
    data.file[data.file=="Medium"] = "DMSO"
	data.values = data.file[,colnames(data.file) %in% model.structure$names]

	# Means of basal activity of the network and of the blank fixation of the antibodies

	unstim.values = colMeans(data.values[data.file$type=="c",])
	blank.values = colMeans(data.values[data.file$type=="blank",])

	# Calculates the mean and standard deviation for each condition 
	mean.values = aggregate(as.list(data.values),by=data.file[,1:3],mean);
	sd.values = aggregate(as.list(data.values),by=data.file[,1:3],sd);
    print(mean.values)

	### CALCULATE ERROR MODEL

	# Define the lower and default error threshold
	min.cv=0.1;     # parameters: minimal cv=0.1
	default.cv=0.3; # parameters: default cv if there are only 2 replicates

	# Calculate error percentage
	cv.values = sd.values[4:dim(sd.values)[2]] / mean.values[4:dim(sd.values)[2]];
	# Values to close to the blank are removed because the error is not due to antibody specific binding
	cv.values[!mean.values[,4:dim(mean.values)[2]] > 2 * matrix(rep(blank.values,each=dim(mean.values)[1]), nrow=dim(mean.values)[1])] = NA;
	
	# Generation of error percentage, one cv per antibody, default.cv if there is only two replicate to calculate the cv
	cv = colMeans(cv.values,na.rm=TRUE)
	cv[cv<min.cv] = min.cv;
    cv[is.nan(cv)|is.na(cv)]=default.cv;
    
    if (FALSE) { #"Multiline comment"
    for (i in 1:dim(cv.values)[2]) {
        count = 0;
        for (j in 1:dim(cv.values)[1]) {
            if (!is.na(cv[i][j]) | !is.nan(cv[i][j])) {
                count = count + 1;
            }
        }
        if (count <= 2) {
            cv[i] = default.cv
            print("Defaulted");
        }
    }
    }

    if (verbose) {
        print("Error model :");
        for (i in 1:length(cv)) {
            print(paste(colnames(data.values)[i], " : ", cv[i]));
        }
    }

# Separate values and perturbation
	data.stim = mean.values[mean.values$type=="t",4:dim(mean.values)[2]];
	data.perturb = mean.values[mean.values$type=="t",1:2]

# Generate the error model
	error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])+matrix(rep(cv,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])*data.stim

### SET UP DATA OBJECT

    data=new(fitmodel::Data);
    data$set_unstim_data (matrix(rep(unstim.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]));
    data$set_scale( data$unstim_data );
    data$set_stim_data( as.matrix(data.stim) );
    data$set_error( as.matrix(error ));


### EXTRACT EXPERIMENTAL DESIGN

# Exctraction of stimulate, inhibited and measured nodes
    stim.nodes=as.character(unique(mean.values$stimulator[mean.values$stimulator %in% model.structure$names]));
    inhib.nodes=as.character(unique(mean.values$inhibitor[mean.values$inhibitor %in% model.structure$names]));
    measured.nodes=colnames(data.stim);

## Identification of nodes with basal activity
    basal.activity=as.character(read.delim(basal_activity,header=FALSE)[,1]);

# Inhibition and stimulation vectors for each experiment
    stimuli=matrix(0,ncol=length(stim.nodes),nrow=dim(data.perturb)[1])
    for (i in 1:length(stim.nodes)) {
	    stimuli[data.perturb$stimulator==stim.nodes[i],i]=1;
    }
    if (verbose) {
        print("Stimulated nodes");
        print(stim.nodes);
        print(stimuli);
    }
    inhibitor=matrix(0,ncol=length(inhib.nodes),nrow=dim(data.perturb)[1])
    for (i in 1:length(inhib.nodes)) {
	    inhibitor[data.perturb$inhibitor==inhib.nodes[i],i]=1;
    }
    if (verbose) {
        print("Inhibited nodes");
        print(inhib.nodes);
        print(inhibitor);
    }

# Experimental design
	expdes=getExperimentalDesign(model.structure,stim.nodes,inhib.nodes,measured.nodes,stimuli,inhibitor,basal.activity);

# Model inputs structure
    model_description = list()
    model_description$design = expdes;
    model_description$structure = model.structure;
    model_description$data = data;

    return(model_description);
}


draw_profiles <- function(model_description=NULL, trace_relation=FALSE)
{
# Finds a minimal model that fits the data and gives its parameters
# Requires the fitmodel package

### MODEL SETUP
    expdes = model_description$design;
    model.structure = model_description$structure;
    data = model_description$data;

    model = new(fitmodel::Model);
    model$setModel(expdes, model.structure);
    if (debug) { # Debug only
        print(model.structure$names);
        print("Adjacency matrix :");
        print(model.structure$adjacencymatrix);
    }

### INITIAL FIT
    print ("Initializing the model parameters…")
    params = c();
    residuals = c();
    nb_samples = 40;
# Random initializations with Latin Hypercube Sample to find a global maximum fit
    samples = qnorm(randomLHS(nb_samples, model$nr_of_parameters()));
    for (i in 1:nb_samples) {
        result = model$fitmodel( data,samples[i,] )
        residuals = c(residuals,result$residuals);
        params = cbind(params,result$parameter)
    }
    if (debug) {
        print(sort(residuals))
    }
    
# Choice of the best fit
    initparams = params[,order(residuals)[1]]
    initresidual = residuals[order(residuals)[1]]
    initial.response = model$getLocalResponseFromParameter( initparams )
    adj = model.structure$adjacencymatrix
    rank = model$modelRank()

    # Plot of the profile likelihood for each path
    init_params = model$getParameterFromLocalResponse(initial.response$local_response, initial.response$inhibitors);
    print(paste(length(init_params), "paths to evaluate"));
    identifiables = model$getParametersLinks();
    #print(identifiables);
    for (path in 1:length(init_params)) {
        lprofile = model$profileLikelihood(data, init_params, path, 1000, 0.01);
        lprofile$residuals[path, lprofile$residuals[path,] >= 2*initresidual] = 2*initresidual; # We do not display very high residual value, as they extend the y-axis which hides lower values

        # Functionnal relations
        pdf(paste("Path_", path, "_profile_likelihood.pdf", sep=""));
        if (trace_relation) {
            par(2, 1);
            title = c(0);
            palette(rainbow(length(init_params)));
            plot(0, 0, ylim=range(lprofile$residuals[-path,]), xlim=range(lprofile$explored)+ c(0, 0.2*(lprofile$explored[length(lprofile$explored)]-lprofile$explored[1])), bty="l");
            for (i in 1:dim(lprofile$residuals)[1]) {
                if (i != path) {
                    lines(lprofile$explored, lprofile$residuals[i,], type="l", col=i);
                    title = c(title, identifiables[i]);
                }
            }
            legend(range(lprofile$explored)[2], range(lprofile$residuals[-path,])[2], title, col=1:length(init_params), lty=1, xpd=T, bty="n");
        }

        # Profile likelihood and identifiability thresholds
        plot(lprofile$explored, lprofile$residuals[path,], type="l");
        lines( lprofile$explored, rep(lprofile$thresholds[1], length(lprofile$explored)), lty=2, col="grey" );
        lines( lprofile$explored, rep(lprofile$thresholds[2], length(lprofile$explored)), lty=2, col="grey" ); # Could be accelerated with two points instead of hundreds
        lines( rep(init_params[path], length(-5:100)), (1 + -5:100/100) * initresidual, col="red");
        title(main = paste("Profile likelihood of", path));
        dev.off()

        print(paste("Path", identifiables[path], "profile likelihood plotted, parameter value =", init_params[path] ));
    }
    
    if (FALSE) { # Verbose ?
        model.structure$setAdjacencymatrix(adj);
        model$setModel(expdes, model.structure);

        print("Parameters :");
        parameters = model$getParameterFromLocalResponse(initial.response$local_response, initial.response$inhibitors); # Identifiables (combination of paths)
        print(parameters);

        print("Response matrix : ");
        print(model.structure$names)
        print(model.structure$adjacencymatrix);

        local_response = model$getLocalResponseFromParameter(model$fitmodel(data, parameters)$parameter);
        print("Local response : ");
        print(local_response$local_response)
        print("Inhibitors :");
        print(local_response$inhibitors);
    }

}



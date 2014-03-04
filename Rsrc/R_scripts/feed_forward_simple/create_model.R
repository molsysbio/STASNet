## Two functions that simplify data parsing and model creation and fitting

# Global variable to have more outputs
verbose = FALSE;
debug = TRUE;

create_model <- function(model.links="links", data.stimulation="data", basal_activity = "basal.dat")
{
# Creates a parametrized model from an experiment file and the network structure
# It requires the file network_reverse_engineering-X.X/r_binding/fitmodel/R/generate_model.R of the fitmodel package
# model.links should be an adjacency list file representing the network
# Experiment file data.stimulation syntax should be as follows, with one line per replicate
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
    if (length(inhib.nodes) > 0) { # Usefull for artificial networks
        for (i in 1:length(inhib.nodes)) {
	        inhibitor[data.perturb$inhibitor==inhib.nodes[i],i]=1;
        }
        if (verbose) {
            print("Inhibited nodes");
            print(inhib.nodes);
            print(inhibitor);
        }
    }

# Experimental design
	expdes=getExperimentalDesign(model.structure,stim.nodes,inhib.nodes,measured.nodes,stimuli,inhibitor,basal.activity);

# Model inputs structure
    model_description = list()
    model_description$design = expdes;
    model_description$structure = model.structure;
    model_description$data = data;
    model_description$fileName = data.stimulation;

    return(model_description);
}


# Finds a minimal model that fits the data and gives its parameters
# Requires the fitmodel package
profile_likelihood <- function(model_description=NULL, trace_relation=FALSE)
{
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
    rank = model$modelRank();

    print("Model simulation :")
    print(model$simulate(data, initparams)$prediction)

    init_params = model$getParameterFromLocalResponse(initial.response$local_response, initial.response$inhibitors);# Useless
    print(paste(length(init_params), "paths to evaluate"));
    identifiables = model$getParametersLinks();
    # Collection of the non identifiables parameters
    ni_profiles = list()
    i_profiles = list();
    for (path in 1:length(init_params)) {
        lprofile = model$profileLikelihood(data, init_params, path, 1000, 0.01);
        lprofile$residuals[path, lprofile$residuals[path,] >= 2*initresidual] = NA;
        print(paste("Parameter", path, "decided"));
        
        # Parameters are non identifiable if their profile likelihood does not reach the low threshold on both sides of the minimum 
        if (lprofile$lowt_upper && lprofile$lowt_lower) {
            i_profiles[[length(i_profiles)+1]] <- lprofile;
        }
        else {
            ni_profiles[[length(ni_profiles)+1]] <- lprofile;
        }
    }

    sorted_profiles = list(i_profiles, ni_profiles);
    print("Profiles extracted and sorted.");

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

    return(sorted_profiles);
}


# Plots the functionnal relation between each non identifiable parameter and their profile likelihood
ni_pf_plot <- function(sorted_profiles, initresidual=0, data_name="default") {
    i_profiles = sorted_profiles[[1]];
    ni_profiles = sorted_profiles[[2]];

# Non identifiables
    nbni = length(ni_profiles);
    print(paste(nbni, "non identifiable paths"));
    pdf(paste("NIplot_", data_name, ".pdf", sep=""), height=2*nbni+1, width=2*nbni+1);
    #limx = i_profiles[[1]]$
    if (nbni > 0) {
        margin = c(2, 2, 0.5, 0.5);
        par(mfcol=c(nbni, nbni), mar=margin, lab=c(3, 3, 4));
        for (ni in 1:nbni) {
            # Set the y margins
            if (ni == 1) {margin[2]=4;}
            else {margin[2]=2;}

            for (j in 1:nbni) {
                # Set the x margins
                if (j == nbni) {margin[1]=4;}
                else {margin[1]=2;}

                par(mar=margin);
                plot(ni_profiles[[ni]]$explored, ni_profiles[[ni]]$residuals[ni_profiles[[j]]$pathid,], xlab="", ylab="", type="l", col=abs(ni-j)+1);
                # Print labels on the right and bottom
                if (ni == 1) { title(ylab=ni_profiles[[j]]$path, las=2); }
                if (j == nbni) { title(xlab=ni_profiles[[ni]]$path); }
                if (j == ni) {
                    # Could be accelerated with two points instead of hundreds
                    lines( ni_profiles[[ni]]$explored, rep(ni_profiles[[ni]]$thresholds[1], length(ni_profiles[[ni]]$explored)), lty=2, col="grey" );
                    lines( ni_profiles[[ni]]$explored, rep(ni_profiles[[ni]]$thresholds[2], length(ni_profiles[[ni]]$explored)), lty=2, col="grey" );
                    if (initresidual != 0) { lines( rep(ni_profiles[[ni]]$value, length(-5:100)), (1 + -5:100/100) * initresidual, col="red"); }
                }
            }
            print(paste("Non identifiable path", ni_profiles[[ni]]$path, "plotted"));
        }
    }

# Identifiables
    nbid = length(i_profiles);
    print(paste(nbid, "identifiable paths"));
    par(mfcol=c(1, 2), mar=rep(5, 4));
    if (nbid > 0) {
        for (id in 1:nbid) {
            plot(i_profiles[[id]]$explored, i_profiles[[id]]$residuals[i_profiles[[id]]$pathid, ], type="l", sub=paste(i_profiles[[id]]$path, "profile"));

            plot(1, type="n", xlim=range(i_profiles[[id]]$explored), ylim=range(i_profiles[[id]]$residuals[-i_profiles[[id]]$pathid], na.rm=T), );
            #print(range(i_profiles[[id]]$residuals[-i_profiles[[id]]$pathid], na.rm=T));
            for (i in 1:dim(i_profiles[[id]]$residuals)[1]) {
                if (i != id) {
                    lines(i_profiles[[id]]$explored, i_profiles[[id]]$residuals[i, ], sub="Functionnal relation");
                }
            }
            print(paste("Identifiable path", i_profiles[[id]]$path, "plotted") );
        }
    }

    dev.off();

}



//SLiM model exploring how pesticide application strategy
//affects the delay of resistance evolution
initialize() {
    initializeSLiMModelType("nonWF");
    defineConstant("Strategy", st); //kind of management: 0 loco (not used),  1 combo, 2 mosaic, 3 rotation, 4 sequential
    defineConstant("RepMode", rm); //reproductive mode: 0 sexual, 1 asexual, 2 alternating
    defineConstant("roTime", rt); //period of periodic strategy
    defineConstant("K", 6000); // total population carrying capacity
    defineConstant("Nsub", 30); // number of subpopulations 
    defineConstant("subK", asInteger(round(K/Nsub))); // subpop carrying capacity
    defineConstant("RefP", rp); //proportion of refugial patches 
    defineConstant("Cost", cos); //cost factor for resitant phenotypes
    defineConstant("e", 0.0); //chance per generation of subpopulation extinction
    defineConstant("m", mi); //migration rate 
    defineConstant("B", 1.0); //mean and variance of Poisson distro of fecundity
    defineConstant("O1", 5.0); // optimum for resistance trait 1 
    defineConstant("O2", 5.0); //optimum for resistance trait 2 
    defineConstant("omega", ss); //weakness of selection
    defineConstant("esigma", 0.2); //variance of environ. effect on phenotypes
    defineConstant("ksigma", 0.2); //adds noise to the optimum for resistance traits
    defineConstant("QTL_mu", c(0, 0)); // pleio effect means
    defineConstant("QTL_cov", c(cv)); // pleio effect covariance
    defineConstant("QTL_sigma", matrix(c(0.5,QTL_cov,QTL_cov,0.5), nrow=2));
    initializeMutationType("m1", 0.5, "n", 0.0, 0.0); // pleiotropic QTLs
    initializeMutationType("m2", 0.5, "n", 0.0, 0.5); // CN t1 QTLs
    initializeMutationType("m3", 0.5, "n", 0.0, 0.5); // CN t1 QTLs
    initializeGenomicElementType("g1", c(m1,m2,m3), c(0.1,1.0,1.0));
    initializeGenomicElement(g1, 0, 99999);
    initializeMutationRate(1e-6);
    initializeRecombinationRate(1e-8);
    m1.convertToSubstitution = F;
    m2.convertToSubstitution = F;
	m3.convertToSubstitution = F;
	m1.mutationStackPolicy = "l";
    m2.mutationStackPolicy = "l";
	m3.mutationStackPolicy = "l";
}

//##########################Mis en place!#################################

//Do cloning (1), or monogamous mating with variable fecundity (0), or alternation between those two (2)
1: reproduction() {
    litter = rpois(1, B);
    if (RepMode==0){
        mate = subpop.sampleIndividuals(1);
        if (mate.size()){ // check to make sure there was a mate
            for (j in seqLen(litter)){
                subpop.addCrossed(individual, mate);
            }
        }
    }
    if (RepMode==2){
        if (sim.generation-1 % 10 == 0){
            mate = subpop.sampleIndividuals(1);
            if (mate.size()){
                for (j in seqLen(litter)){
                    subpop.addCrossed(individual, mate);
                }
            }
        } else {
            for (j in seqLen(litter)){
                subpop.addCloned(individual);
            }
        }
    }
    else {
        for (j in seqLen(litter)){
            subpop.addCloned(individual);
            }
    }
}

//This is a quantitative genetic model,
//so turn off direct allele fitness effects
fitness(m1) { return 1.0; }
fitness(m2) { return 1.0; }
fitness(m3) { return 1.0; }

//since we are allowing for pleiotropies, we make a special
//mutation generator that samples t1 and t2 effects from
//a bivariate normal distribution
mutation(m1){
	// draw mutation effect for new m2 pleiotropic mutation
	effects = rmvnorm(1, QTL_mu, QTL_sigma);// * Mutation_var; this to scale 
    mut.setValue("e1", effects[0]);
    mut.setValue("e2", effects[1]);
	return T;
}

//Add a Npop subpopulations
//We also creat a variable to -- when using the sequential strategy -- register
//if the first treatment has stopped working
//and a variable to keep track of which poison we're using in a rotational strategy
1 early() {
    for (i in 1:Nsub){
        sim.addSubpop(i, subK);
    }
    sim.setValue("SmokeEm", 0);
    sim.setValue("HotPoison", 0);
    sim.setValue("baddies", 0);
}

//do non-overlapping generations
1: early() {
    inds = sim.subpopulations.individuals;
    inds[inds.age > 0].fitnessScaling = 0.0;
}

//set the control strategy
1 early() {
    //some patches are untreated refuges
    swubs = seq(0,Nsub-1,by=1);
    refuges = sample(swubs, asInteger(round(RefP*Nsub)));
    //kinda ugly, but this gives us a vector of treated field indices
    blast = rep(T, Nsub);
    for (i in refuges){
        blast[i] = F;
    }
    treated = swubs[blast];
    Ntreated = size(treated);
    
    //Set which fields get which treatment, depending on Strategy
    unos = rep(0, Ntreated); //vector specifying which patches will get treatment 1
    doses = rep(0, Ntreated); //vector of which patches that get treatment 2
    if (Strategy == 1){ //use combo approach
        unos = rep(1, Ntreated);
        doses = rep(1, Ntreated);
    } if (Strategy == 2){ //use the mosiac approach
        unos = rbinom(Ntreated, 1, 0.5);
    } else { //use either the rotation or sequential approach
        unos = rep(1, Ntreated);
    }
    //plop those into global variables
    sim.setValue("unos", unos);
    sim.setValue("doses", doses);
    sim.setValue("refuges", refuges);
    sim.setValue("Ntreated", Ntreated);
    sim.setValue("treated", treated);
    sim.setValue("mW", 1.0);
    sim.setValue("rN", 800);
}

//In each generation we'll just have m*N individuals move randomly from 
//their current patch to any other. We also allow for random extinction
//of subpopulaitons. This follows the recipe in section 16.5 of the SLiM Manual
1: early(){
    nIndividuals = sum(sim.subpopulations.individualCount);
    if (nIndividuals < 1){
        catn("Pest Extinction!!");
        sim.simulationFinished();
    }
    nMigrants = rpois(1, nIndividuals * m);
    migrants = sample(sim.subpopulations.individuals, nMigrants);
    for (migrant in migrants){
        do dest = sample(sim.subpopulations, 1);
        while (dest == migrant.subpopulation);
        dest.takeMigrants(migrant);
    }
    // density-dependence and random extinctions
    for (subpop in sim.subpopulations){
        if (runif(1) < e){
            subpop.fitnessScaling = 0.0;
        } else {
            subpop.fitnessScaling = min(subK / subpop.individualCount, 1.5);
        }
    }
}

//##########################Now Get Into It!#################################

//After 100 gens, we start applying control strategies
50: early() { 
    optimum1 = O1 + rnorm(1, 0, ksigma);
    optimum2 = O2 + rnorm(1, 0 , ksigma);
    //retrieve some variables set in generation 1.
    unos = sim.getValue("unos");
    doses = sim.getValue("doses");
    refuges = sim.getValue("refuges");
    Ntreated = sim.getValue("Ntreated");
    treated = sim.getValue("treated");
    //in refuges there's only costs for resistance alleles
    for (i in refuges){
        subby = sim.subpopulations[i];
        inds = subby.individuals;
        for (ind in inds){
            //ind = inds[i-1];
            muts1 = ind.genomes.mutationsOfType(m1);
            muts2 = ind.genomes.mutationsOfType(m2);
            muts3 = ind.genomes.mutationsOfType(m3);
            t1genotype = size(muts1) ? sum(muts1.getValue("e1")) else 0.0;
            c1genotype = size(muts2) ? ind.sumOfMutationsOfType(m2) else 0.0;
            t1phenotype = t1genotype + c1genotype + rnorm(size(t1genotype),0,esigma);
            t2genotype = size(muts1) ? sum(muts1.getValue("e2")) else 0.0;
            c2genotype = size(muts3) ? ind.sumOfMutationsOfType(m3) else 0.0;
            t2phenotype = t2genotype + c2genotype + rnorm(size(t2genotype),0,esigma);
            Costs = abs(Cost*t1phenotype) + abs(Cost*t2phenotype);
            Cs = min(1.0, Costs);
            ind.fitnessScaling = 1 - Cs;
        }
    }
    //if the strategy is loco then we update the treated fields each gen
    if (Strategy == 0){ //use loco approach
        unos = rbinom(Ntreated, 1, 0.5);
        doses = rbinom(Ntreated, 1, 0.5);
    } 
    if (Ntreated > 0){
        for (i in 1:Ntreated){
            j = treated[i-1];
            uno = unos[i-1];
            dos = doses[i-1];
            inds = sim.subpopulations[j].individuals;
            t1adaptation = rep(1.0,size(inds));
            t2adaptation = rep(1.0,size(inds));
            counter = 0;
            for (ind in inds){
                //ind = inds[i-1];
                muts1 = ind.genomes.mutationsOfType(m1);
                muts2 = ind.genomes.mutationsOfType(m2);
                muts3 = ind.genomes.mutationsOfType(m3);
                t1genotype = size(muts1) ? sum(muts1.getValue("e1")) else 0.0;
                c1genotype = size(muts2) ? ind.sumOfMutationsOfType(m2) else 0.0;
                t1phenotype = t1genotype + c1genotype + rnorm(size(t1genotype),0,esigma);
                t1deviation = optimum1 - t1phenotype;
                fitnessFunctionMax = dnorm(0.0, 0.0, omega);
                t1adapt = dnorm(t1deviation, 0.0, omega) / fitnessFunctionMax;
                t1adaptation[counter] = t1adapt;
                
                t2genotype = size(muts1) ? sum(muts1.getValue("e2")) else 0.0;
                c2genotype = size(muts3) ? ind.sumOfMutationsOfType(m3) else 0.0;
                t2phenotype = t2genotype + c2genotype + rnorm(size(t2genotype),0,esigma);
                t2deviation = optimum2 - t2phenotype;
                fitnessFunctionMax = dnorm(0.0, 0.0, omega);
                t2adapt = dnorm(t2deviation, 0.0, omega) / fitnessFunctionMax;
                t2adaptation[counter] = t2adapt;
                counter = counter + 1;
                Costs = abs(Cost*t1phenotype) + abs(Cost*t2phenotype);
                Cs = min(1.0, Costs);
                ind.fitnessScaling = 1 - Cs;
            }
            
            if (Strategy == 0 | Strategy == 1){ //do loco or combo
                if (uno==1 & dos==1){
                    inds.fitnessScaling = (t2adaptation * t1adaptation);
                } if (uno==1 & dos==0){
                    inds.fitnessScaling = t1adaptation;
                } if (uno==0 & dos==1){
                    inds.fitnessScaling = t2adaptation;
                } else {
                    //no control in this patch this time
                    inds.fitnessScaling = inds.fitnessScaling;
                }
            } if (Strategy == 2){ //do mosaic
                if (uno == 1){
                    inds.fitnessScaling = t1adaptation;
                } else {
                    inds.fitnessScaling = t2adaptation;
                }
            } if (Strategy == 3){ //do rotation
                //an implementation in which we rotate at a specified number of generations
                hot = sim.getValue("HotPoison");
                if ((sim.generation-1) % roTime == 0){ //switch it
                    hot = abs(hot - 1);
                    sim.setValue("HotPoison", hot);
                }
                if (hot == 0){
                    inds.fitnessScaling = t1adaptation;
                } else {
                    inds.fitnessScaling = t2adaptation;
                }
                
            } if (Strategy == 4){ //do sequential
            
                if (sim.getValue("SmokeEm") == 0){
                    inds.fitnessScaling = t1adaptation;
                } else {
                    inds.fitnessScaling = t2adaptation;
                }
            }

        }
     }
}

//Just check on mean fitnesses
50: early(){
    if ((sim.generation-1)%10==0){
        W = sim.getValue("mW");
        treats = sim.getValue("treated");
        for (chi in treats){
            pespo = sim.subpopulations[chi];
            mw = mean(pespo.individuals.fitnessScaling);
            W = c(W, mw);
            sim.setValue("mW", W);
        }
        foo = size(sim.subpopulations.individuals);
        sim.setValue("rN", foo);
    }
}

//OK. Now keep track of some stuff, and throw in the towel when the 
//pest population has reached > 25% of the carrying capacity of treated fields.
50: late(){
    if ((sim.generation -1) % 10 == 0){
        //determine the threashold pest pop
        subKi = K/Nsub;
        nTreated = sim.getValue("Ntreated");
        thresh = 0.25 * (nTreated*subKi);
        
        inds = sim.subpopulations.individuals;
        for (ind in inds){
            muts1 = ind.genomes.mutationsOfType(m1);
            muts2 = ind.genomes.mutationsOfType(m2);
            muts3 = ind.genomes.mutationsOfType(m3);
            t1genotype = size(muts1) ? sum(muts1.getValue("e1")) else 0.0;
            c1genotype = size(muts2) ? ind.sumOfMutationsOfType(m2) else 0.0;
            t1tot = t1genotype + c1genotype;
            t2genotype = size(muts1) ? sum(muts1.getValue("e2")) else 0.0;
            c2genotype = size(muts3) ? ind.sumOfMutationsOfType(m3) else 0.0;
            t2tot = t2genotype + c2genotype;
            ind.setValue("t1g", t1tot);
            ind.setValue("t2g", t2tot);
        }
        t1s = mean(inds.getValue("t1g"));
        t2s = mean(inds.getValue("t2g"));
        mu1 = mean(t1s);
        mu2 = mean(t2s);
        pests = 0;
        mw = mean(sim.getValue("mW"));
        for (ips in sim.getValue("treated")){
            farm = sim.subpopulations[ips];
            pests = pests + farm.individualCount;
        }
        catn(c(sim.generation-1,mu1,mu2,sum(sim.subpopulations.individualCount),pests,mw,sim.getValue("rN")));
        bads = sim.getValue("baddies");
        sim.setValue("baddies", c(bads, sim.subpopulations.individualCount));
        //look at top resistance phenos
        delta1 = min(abs(O1-t1s));
        delta2 = min(abs(O2-t2s));
        closeEnough = 0.5*abs(O1);
        //if (pests > thresh){
        if (delta1 < closeEnough | delta2 < closeEnough ){
            catn("Control failure!!");
            city = mean(sim.getValue("baddies"));
            sim.simulationFinished();
            writeFile("burnOut.csv", paste(c(Strategy,QTL_cov,RepMode,roTime,m,omega,Cost,city,(sim.generation-1)),sep=','), append=T);
        } 
    }
}

50000 late() { 
    city = mean(sim.getValue("baddies"));
    sim.simulationFinished(); 
    writeFile("burnOut.csv", paste(c(Strategy,QTL_cov,RepMode,roTime,m,omega,Cost,city,(sim.generation-1)),sep=','), append=T);
}

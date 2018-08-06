#############################
#----- Constants -------#
#############################

g0=9.81 #gravity
cp=1005 #specific heat of air 
lv=2.5*10^6 #latent heat of vaporization
p0=100000 #reference pressure
Rd=287.04 #gas constant
Rm=461.5 
ep2=Rm/Rd-1
th00=289 #reference temperature
tmelt=273.16 #melting point of water


#############################
#----- Graphic details -----#
#############################

#letters
lettering=paste0("(",letters,")")

#octant names and colors
octname=c("Entrainment shell","Mixing updraft","Ascending shell","Updraft","Entrainment","Mixing shell","Downdraft","Subsiding shell")
PALOCT=c("lightpink","yellow","orange2","red3","palegreen","darkgreen","navy","royalblue")


#############################
#----- Spatial tools -------#
#############################

# Vertical mean (slab-average)
spatial.mean<-function(field) {
	outfield=colMeans(field,dims=2,na.rm=T)
        return(outfield)
}

# Slab average anomalies
eddy.mean<-function(field) {
	outfield=sweep(field,c(3),spatial.mean(field),"-")
	return(outfield)
}

#detect level
whicher<-function(axis,number) {
	        out=which.min(abs(axis-number))
        return(out)
}

#############################
#----- Thermodynamics ------#
#############################

# saturation pressure of vapour at defined temperature
esl<-function(t) {
	c0=0.6105851e+03;
	c1=0.4440316e+02;
	c2=0.1430341e+01;
	c3=0.2641412e-01;
	c4=0.2995057e-03;
	c5=0.2031998e-05;
	c6=0.6936113e-08; 
	c7=0.2564861e-11;
	c8=-.3704404e-13;
	x=t-tmelt
    esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
    return(esl)
}

# Saturation water vapour mixing ratio for specific temperature and pressure
qsat<-function(t,p) {
	qsat=0.622*esl(t)/(p-esl(t))
	return(qsat)
}

# Exner function
exner<-function(p,p0=10^5) {
	exner=(p/p0)^(Rd/cp)
	return(exner)
}

# Saturation adjustment from liquid potential temperature, 
# total water vapour mixing ratio and pressure
# it uses optimization function by R
full.satadj<-function(theta_l,q_tot,p,lv=2.5*10^6) {

        # cost function to be minimized: for a give temperature
        # the difference between the liquid water mixing ratio 
        # given by liquid potential temperature equation and 
        # water vapour saturation should be smaller as possible
        core.satadj<-function(tx) {
                theta=tx/exner(p)
                l1=(theta-theta_l)*(cp/lv)
                l2=q_tot-qsat(tx,p)
                delta=abs(l1-l2)
                return(delta)
        }

        # minimized the cost function
        conv=optimize(core.satadj,lower=200,upper=theta_l)

	#build userful variables
        tx=conv$minimum
        q_sat=qsat(tx,p)
        theta=tx/exner(p)
        l1=q_tot-q_sat
        l2=(theta-theta_l)*(cp/lv)
        l=(l1+l2)/2
	# if the obtained liquid water is negative a few adjustment
	# should be made
        if (l<0) {
        	l=0
        	theta=theta_l
        	tx=theta*exner(p)
        	q_sat=qsat(tx,p)
        }
        
	#virtual potential temperature      
        theta_v=theta*(1+0.61*(q_tot-l)-l)
	buoyancy=(theta_v-th00)/th00*g0


        #return values for all thermodynamic variables
        out=list(theta_l=theta_l,theta_v=theta_v,theta=theta,
        tx=tx,q_tot=q_tot,q_v=q_tot-l,q_sat=q_sat,l=l,p=p,buoyancy=buoyancy)
        return(out)

}

mixing.curve<-function(t1,q1,p1,t2,q2,p2,LV=2.5*10^6) {
	
	#define k parameter for CTEI
	k=1+(t2-t1)/(LV/cp*(q2-q1))

	#define mixing function
	mixing<-function(alpha,lv=LV) {
	
		tm=alpha*t1+(1-alpha)*t2
		qm=alpha*q1+(1-alpha)*q2
		pm=alpha*p1+(1-alpha)*p2

		sat=full.satadj(tm,qm,pm,lv=LV)
		return(sat$theta_v)
	}

	#find minimum buoyancy reversal
	ki=optimize(mixing,lower=0,upper=1)
	mixing_curve=c(mixing(0),mixing(ki$minimum),mixing(1))
	ratio=c(0,round(ki$minimum,4),1)

	#compute buoyancy
	buoyancy=(mixing_curve-th00)/th00*g0
	theta_l=c(t2,(t1-t2)*ratio[2]+t2,t1)
	q_tot=c(q2,(q1-q2)*ratio[2]+q2,q1)
	out=list(ratio=ratio,buoyancy=buoyancy,k=k,ki_s=round(ki$minimum,4),
	theta_v=mixing_curve,theta_l=theta_l,q_tot=q_tot)
	return(out)
}


#############################
#------ NetCDF tools -------#
#############################

# lasthours.file
# from a netcdf file extract time array and identify last hour records
lasthour.file<-function(filename,timeblock=3600) {
	require(ncdf4)
	filein=nc_open(filename)
	time=ncvar_get(filein,"time")
	lasthour_start=length(time)-timeblock/diff(time)[1]+1
	lasthour_end=lasthour_start+timeblock/diff(time)[1]-1
	return(list(time=time,lasthour=lasthour_start:lasthour_end))
	nc_close(filein)
}

#############################
#------ Derivatives --------#
#############################


der.zeta<-function(ics,ipsilon,vertical,field,jump=2) {
	
	dims=dim(field)
	if (is.null(dims)) {d=1}
	if (!is.null(dims)) {d=length(dims)}
	outfield=field*NA

	if (d==1) {
		delta=diff(field,jump)
		increment=(diff(vertical,jump))
		if (jump%%2==0) {outfield[(jump/2+1):(length(vertical)-(jump/2))]=delta/increment}
		if (jump%%2==1) {outfield[(1+jump%/%2):(length(vertical)-ceiling(jump/2))]=delta/increment}
	}


	if (d==3) {
		for (i in 1:length(ics))
		{for (j in 1:length(ipsilon))
		{
			#differenzial of the field
			delta=diff(field[i,j,],jump)

			#compute the increment
			increment=(diff(vertical,jump))

			#compute the final matrix, considering the differenting window adopted
			if (jump%%2==0) {outfield[i,j,(jump/2+1):(length(vertical)-(jump/2))]=delta/increment}
			if (jump%%2==1) {outfield[i,j,(1+jump%/%2):(length(vertical)-ceiling(jump/2))]=delta/increment}
		}
		}
	}

	if (d==2) {
		for (j in 1:length(ipsilon))
		{
			#differenzial of the field
			delta=diff(field[j,],jump)

			#compute the increment
			increment=(diff(vertical,jump))

			#compute the final matrix, considering the differenting window adopted
			if (jump%%2==0)
			{outfield[j,(jump/2+1):(length(vertical)-(jump/2))]=delta/increment}
			if (jump%%2==1)
			{outfield[j,(1+jump%/%2):(length(vertical)-ceiling(jump/2))]=delta/increment}
		}

	}

	return(outfield)
}



#############################
#--- Coherent Structures ---#
#############################

# BruteForce.fast function
# estimates the vertical integral of heat flux for downdraft and entrainment octants
bruteforce.fast<-function(scft_conv=0) {

	#start
	#t0=proc.time()

	# define arrays
	octo=w_prime*0
	scft_star=scft-scft_conv

	#timing 1
	#t1b=proc.time()
	#printv(t1b-t0)

	#use this predeclaration to save a few time	
	pp<-(w_prime<0 & scbl_prime>0)
	octo[pp & scft_star<0]=7
	octo[pp & scft_star>0]=5

	#timing 2
	#t1=proc.time()
	#print(t1-t1b)

	# OLD METHOD
	#counting cases
	#octperc5=apply(octo,3,function(x) (length(x[x==5])))/np
	#octperc7=apply(octo,3,function(x) (length(x[x==7])))/np

	# timing 3
	#t2=proc.time()
	#print(t2-t1)

	#heatflux
	#pdown=pentrain=2:length(lev)*NA
	#for (z in 1:(length(lev)-1)) {
	#	heatflux0=heatflux[,,z]
	#	if (octperc5[z]>0) {pentrain[z]=octperc5[z]*mean(heatflux0[octo[,,z]==5])}
	#	if (octperc7[z]>0) {pdown[z]=octperc7[z]*mean(heatflux0[octo[,,z]==7])}
	#}
	
  	# timing 3
        #t2=proc.time()
        #print(t2-t1)

	# NEW METHOD
        pdown=pentrain=2:length(lev)*NA
        for (z in 1:(length(lev)-1)) {
                heatflux0=heatflux[,,z]
		p5 <- octo[,,z]==5
		p7 <- octo[,,z]==7
		octperc5=length(p5[p5==TRUE])/np
		octperc7=length(p7[p7==TRUE])/np
                if (octperc5>0) {pentrain[z]=octperc5*mean(heatflux0[p5])}
                if (octperc7>0) {pdown[z]=octperc7*mean(heatflux0[p7])}
        }


	# timing 4
	#t3=proc.time()
	#print(t3-t2)

	#final
	sumpdown=sum(pdown*diff(lev),na.rm=T); #print(sumpdown)
	sumpentrain=sum(pentrain*diff(lev),na.rm=T); #print(sumpentrain)
	outfield=sumpdown-sumpentrain

	#assign("pdown",pdown, envir = .GlobalEnv)
	#assign("pentrain",pentrain, envir = .GlobalEnv)
	return(outfield)
}


# Convergence function
# Estimates the values for maximum difference in the heat flux 
convergence<-function(min.digits=3,max.digits=6,method="serial",ncluster=4,savememory=F,verbose=T) {

# check for method availability
# NB: parallelization is not scaling as may be expected
# Zero-order benchmark suggest that with 128x128x210 3D domain on our local Unix cluster
# 1. 54 seconds for serial
# 2. 48 seconds with 4 cores and parSapply
# 3. 40 seconds with 2 cores and foreach
# 4. 30 seconds with 4 cores and foreach
# 5. 32 seconds with 8 cores and foreach
# so that foreach is the suggested method if parallelization is available, but using many cores is pointless. 
# Savememory removes a bit of the efficiency but it should be useful with large datasets
if (!any(method==c("serial","foreach","parSapply"))) {
	stop("Only 3 method supported so far: serial, foreach and parSapply")
}

# useful function for verbosity
printv<-function(value) {if (verbose) {print(value)} }

# Verbose information on the convergence
print("Convergence algorithm to maximize delta fluxes")
print(paste("Precision is set from",min.digits,"to",max.digits)) 
print(paste("Convergence mode is",method))
if (method!="serial") {
	print("Parallel mode is on")
	print(paste("Running on",ncluster,"cores"))
	print(paste("Save-memory is",savememory))
}

# First guess, done with optimize function and min.digits precision
print("First guess...")
t0=proc.time()
convergence1=optimize(bruteforce.fast,lower=0,upper=.1,maximum=TRUE,tol=10^-(min.digits))
maximum=convergence1$maximum; objective=convergence1$objective
print(paste("Threshold:",round(maximum,min.digits)," - Maximum: ",round(objective,min.digits)))
printv(proc.time()-t0)
printv("---------------")

# Loop on following increasing precision: the function is unimodal but since it is
# this is done in order to avoid numerical issues in the neighbourhood of the maximum
# the increasing precision is used to reduce the number of intervals required
for (d in min.digits:max.digits) {
        print(paste0("Decimal Loop #",d))
	#printv(mem_used())

	#Number of intervals on which we will run convergence
        span=d+4
        if (d==min.digits) {objective=convergence1$objective; maximum=convergence1$maximum}
        delta=10^(-d)

	# serial convergence with sapply on the -span:span interval
	if (method=="serial") {
                convergence2=sapply(-span:span,function(i) optimize(bruteforce.fast,
				   lower=maximum+(delta*i)-delta/2,upper=maximum+(delta*(i+1))-delta/2,
				   maximum=TRUE,tol=10^(-i)))
                }

	# if using a parallel method, activate the required stuff
        if (method!="serial") {
                require(parallel); require(doParallel)

		# savememory option, for large dataset avoid complete duplication of environments
                if (savememory) {
                        cl <- makeCluster(ncluster)
                        printv("Savememory")
                        clusterExport(cl, varlist=list("heatflux","scbl_prime","w_prime","scft","np","lev","bruteforce.fast"))
                        clusterExport(cl, varlist=list("delta","maximum"),envir=environment())
                } else {
                	cl <- makeCluster(ncluster,type="FORK")
		}

		#activate cluster
                registerDoParallel(cl)
                printv("Parallel convergence...")

		# for each method
		if (method=="foreach") {
			require(foreach)
			
			#define for each function
                	f1<-function(span,maximum,delta) {
                        foreach(i=-span:span,.combine=cbind) %dopar% {
                                p=optimize(bruteforce.fast,lower=maximum+(delta*i)-delta/2,upper=maximum+(delta*(i+1))-delta/2,maximum=TRUE,tol=10^(-i))
                                return(unlist(p))

                        	}
                	}
			
			# run convergence
                	convergence2=f1(span,maximum,delta)
		}

		#sapply method
		if (method=="parSapply") {
                	convergence2=parSapply(cl,-span:span,function(i) optimize(bruteforce.fast,
					       lower=maximum+(delta*i)-delta/2,upper=maximum+(delta*(i+1))-delta/2,
					       maximum=TRUE,tol=10^(-i)))
		}

		#stop cluster for parallel computation
		stopCluster(cl)
        }

	# rearrange outputs
        matrice=matrix(unlist(convergence2),2)
        newobjective=max(c(matrice[2,]))

	# if the new objective is better than before, redfined all the story
        if (newobjective>objective) {
                interval=which.max(c(matrice[2,]))-span-1
                printv(paste0("Found new Maximum at interval ",interval))
                objective=newobjective
                maximum=matrice[1,matrice[2,]==objective[1]]
        } else { 
		printv("Keep old Maximum")
	}
        print(paste("Threshold:",round(maximum[1],d)," - Maximum: ",round(objective[1],d)))
        printv("---------------")
}

# returing results
printv(proc.time()-t0)
print(paste("Running time is",round((proc.time()-t0)[3],2),"seconds"))
return(list(maximum=maximum,objective=objective))
}

# Shape Factor function
shape.factor.fast<-function(lon,lat,field,octvalue,min.area=50000)
{
	require(SDMTools) # for connected component labeling
        #properties as function of the octant
        if (octvalue==5) {
		field[field==5 | field==7]=10; field[field!=10]=0; field[field!=0]=1
		} else {
		field[field!=octvalue]=0; field[field!=0]=1
	}
        field[is.na(field)]=0

        #connected componenent labeling and stima delle aree
        ppp=ConnCompLabel(field)
        #print(quantile(ppp,na.rm=T))
        ncounts=table(factor(ppp,levels=1:max(ppp)))
        #convert area in grid points
        dx=diff(lon)[1]; dy=diff(lat)[1]
        min.area=min.area/(dx*dy)

        #take only elements big enough
        elements=which(ncounts>=min.area & as.numeric(rownames(ncounts))>0)
        if (length(elements)>0) {
        #print(length(elements))
        #print(ncounts[elements])
        	circularity=1:length(elements)*NA
        	for (element in elements) {
                	#isolate element and estimate area
                	a=ppp; a[a!=element]=0
                	area=length(a[a!=0])*dx*dy
                	#extend boundaries
                	att1=lon*0; att2=c(0,lat*0,0)
                	b=cbind(att1,a,att1); b=rbind(att2,b,att2)
                	extlon=c(lon[1]-dx,lon,lon[length(lon)]+dx)
                	extlat=c(lat[1]-dy,lat,lat[length(lat)]+dy)
                	#compute contour lines for perimeter
                	ll=contourLines(extlon,extlat,b,level=element)
                	#identify larger contour line
                	rr=NULL
                	#for (k in 1:length(ll)) {rr=c(rr,length(ll[[k]]$x))}
                	rr=sapply(1:length(ll),function(k) length(ll[[k]]$x))
                	linea=ll[[which.max(rr)]]
                	#compute perimeter as sum of successive distance of points
                	perimeter=0
                	#t1=proc.time()
                	#for (i in 2:length(linea$x))
                	#        {perimeter=perimeter+dist(cbind(linea$x[(i-1):i],linea$y[(i-1):i]))}
                	pre=sapply(2:length(linea$x),function(i) {dist(cbind(linea$x[(i-1):i],linea$y[(i-1):i]))})
                	perimeter=sum(pre)
                	#print(proc.time()-t1)
	
        	        #element circularity
                	circularity[which(element==elements)]=(4*pi*area)/(perimeter^2)
                }
        #weighted circularity
        circ=weighted.mean(circularity,ncounts[elements])
        } else {
        	circ=NA
        }
        return(circ)
}

#compute fluxes contribution
profile.mean<-function(var) {
	if (var=="octo") {
		profile=lev*0
	} else {
		profile=colMeans(get(var),dims=2)
	}
return(profile)
}

#compute octants profiles
profile.octants.new<-function(var,octo,counts) {
	var0=get(var)
	if (var=="massflux" | var=="wtheta" | var=="wqtot") {
		thr=0
	} else {
		thr=0.02
	}
	np=length(var0[,1,1])*length(var0[1,,1])
	if (var=="octo") {
		profile=t(counts)/np*100
	} else {
		profile=sapply(1:length(lev), function(i) {
                	v0=var0[,,i];
                	mV=sapply(1:length(counts[,1]), function(j) { if (counts[j,i]/np>thr) { mean(v0[octo[,,i] == j], na.rm = TRUE) } else {NA} })
                	if (var=="massflux" | var=="wtheta" | var=="wqtot") {counts[,i]/np*mV} else {mV}
        		} )
	profile<-matrix(profile,length(lev),length(counts[,1]),byrow = TRUE)
	}
return(profile)
}



# diagnostic

#memory usage
sort( sapply(ls(),function(x){format(object.size(get(x)),unit="auto")}))


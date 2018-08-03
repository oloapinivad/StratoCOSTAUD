#program for octants and figures for dycoms2_rf01 stratocumulus convection

library(SDMTools) # for connected component labeling
CODEDIR="/home/paolo/StratoCostaud"
source(file.path(CODEDIR,"les_routines.R"))

args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	expcode="tecs00"
	ntimes=seq(1)
} else {
	expcode=args[1]
	ntimes=as.numeric(args[2])
}

#setup flags
do_loading=T
do_computation=T
do_parallel=T
conv_method="foreach"
cores=4
do_savememory=T
do_sink=F

# CONFIGURATION 
model="uclales2"
exptype="dycoms2_rf01"
DIRBASE=file.path("/work/users/paolo",model)
FIGBASE="/home/paolo/scratch/coherent_test"

#vars to be loaded
varlist=c("u","v","w","scbl","scft","qtot","tl","ql","rflx","prss")

# variable names
u_name="u"
w_name="w"
v_name="v"
scbl_name="scbl"
scft_name="scft"
qtot_name="q"
tl_name="tl"
ql_name="l"
prss_name="pr"
rflx_name="rf"

# experiment-dependent properties
if (exptype=="dycoms2_rf01") {
	u0=7; v0=-5.5
}

# extra manipulation
expname=paste(exptype,expcode,sep="_")

# ------------------------------- #
# ------- START PROGRAM --------- # 
# ------------------------------- #

print(paste("Expname:",expname))
for (nt in ntimes) {
	padnt=sprintf("%02d",nt)
	print(paste("Analysing",padnt,"timestep for",expcode))

	# sink: check the goal
	if (do_sink) {
		sinkfile=paste0("/home/paolo/scratch/coherent_test/test_sink_",expcode,"_",padnt,do_parallel,do_savememory,".txt")
		sink(file=sinkfile)
	}

	#dirs
	FIGDIR=file.path(FIGBASE,expname,"evolution",padnt)
	dir.create(FIGDIR,recursive=T)

	if (do_loading) {

	#set ts and ps file
	fileps=paste0(file.path(DIRBASE,expname),"/",expname,".py.ps.nc")
	filets=paste0(file.path(DIRBASE,expname),"/",expname,".py.ts.nc")
	file3d=paste0(file.path(DIRBASE,expname),"/",expname,".py.3d.nc")

	# load density
	filein=nc_open(fileps)
	dn0=ncvar_get(filein,"dn0")
	nc_close(filein)

	# load time series (of clouds)
	varts=c("zcmn","zbmn")
	filein=nc_open(filets)
	for (var in varts) {
		field=ncvar_get(filein,var)
		assign(paste(var),field[length(field)])
	}

	# load 3d field
	filein=nc_open(file3d)
	varnames=paste(varlist,"name",sep="_")
	for (var in varnames) {
        	field=ncvar_get(filein,get(var),start=c(1,1,1,nt),count=c(-1,-1,-1,1))
		#print(str(field))
        	assign(paste0(varlist[which(var==varnames)]),field)
		print(paste("Loading",varlist[which(var==varnames)],"which in is model variable",get(var)))
		rm(field)
	}
        
	lon=ncvar_get(filein,"xt")
       	lat=ncvar_get(filein,"yt")
	lev=ncvar_get(filein,"zt")
	nc_close(filein)

	np=length(lon)*length(lat)
	
	u=u0+u; v=v0+v
	#print(mem_used())
	}

	if (do_computation) {
	print("start computation...")

	#fix water vapour and define companion variables
	th=tl+lv/cp*ql
	tk=th*1/((p0/prss)^(Rd/cp))
	tv=th*(1+0.61*(qtot-ql)-ql)
	#mse=sweep(cp*tk+(qtot-ql)*lv,c(3),g0*lev,"+")

	#qrad
	#qrad=(der.zeta(1,1,lev,rflx[1,1,]))/(dn0[1:length(rflx[1,1,])]*cp)*86400
	qrad3d=-sweep(der.zeta(lon,lat,lev,rflx)*3600,c(3),(dn0[1:length(rflx[1,1,])]*cp),"/")
	#maxrad=max(qrad,na.rm=T) #K/day
	#hrad=lev[which(qrad==max(qrad,na.rm=T))] #node of the hyperbolic tanget: max cooling
	
	#PBL height defition
	zi=lev[which.max(diff(apply(tl,c(3),mean,na.rm=T)))] #end of the inversion layer, maximum gradient of T

	#2D-PBL height definition
	#print("2d PBL definition")
	slevel=0.1
	zi2d=apply(tl,c(1,2),function(x) lev[which.max(diff(x))])
	szi2d=apply(scft,c(1,2),function(x) lev[which(x>slevel)[1]])

	#supplementary definition of PBL height
	#varqtot=apply(qtot_prime,3,sd)^2
	#zilow=lev[whicher(varqtot[50:which.max(varqtot)],0.05*max(varqtot))+50-1]
	#zihigh=lev[whicher(varqtot[which.max(varqtot):length(varqtot)],0.05*max(varqtot))+which.max(varqtot)-1]

	#primes
	print("eddies...")
	for (var in c("w","scbl","scft","tl","qtot")) {
		assign(paste0(var,"_prime"),eddy.mean(get(var)))
	}

	#buoyancy and fluxes and tke
	print("fluxes and tke...")
	massflux=sweep(w_prime,c(3), dn0, "*")
	wtheta=heatflux=w_prime*tl_prime
	wqtot=qtot_prime*w_prime
	buoyancy=eddy.mean((tv-th00)/th00*g0) #model definition
	tke=1/2*(eddy.mean(u)^2+eddy.mean(v)^2+w_prime^2)

	#select boundary layer points
	zz=1:whicher(lev,zi)

	#cleaning useless var
	#print(mem_used())
	print("cleaning memory...")
	rm(u,v,prss,tk,th)
	#print(mem_used())

	#brute force
	print("brute force")

	#define filenames for saved files
	#scft_conv=0
	save_file=file.path(FIGDIR,paste(padnt,expcode,"bruteforce_scft_conv",sep="_"))
	save_profiles=file.path(FIGDIR,paste(padnt,expcode,"profiles",sep="_"))

	# the big guy: convergence algorithm here
	print("Convergence optimization for the Free Tropospheric Scalar...")
	solved=convergence(method=conv_method,ncluster=cores,savememory=do_savememory,verbose=F)
	maxflux=solved$objective
	scft_conv=solved$maximum
	#print(mem_used())
	print("Convergence reached: finalizing...")

	# save the useful files for convergence
	scft_zi=mean(scft[,,zz])
	save(file=save_file,lev,scft_conv,scft_zi,szi2d,zi2d)

	# define octants
	print("octants...")
	octo=w_prime*NA
	scft_star=scft-scft_conv
	octo[w_prime>0 & scft_star>0 & scbl_prime<0]=2
	octo[w_prime>0 & scft_star<0 & scbl_prime>0]=3
	octo[w_prime<0 & scft_star>0 & scbl_prime<0]=6
	octo[w_prime<0 & scft_star<0 & scbl_prime>0]=7
	octo[w_prime>0 & scft_star>0 & scbl_prime>0]=1
	octo[w_prime>0 & scft_star<0 & scbl_prime<0]=4
	octo[w_prime<0 & scft_star>0 & scbl_prime>0]=5
	octo[w_prime<0 & scft_star<0 & scbl_prime<0]=8

	# clean octants
	octo[scft>0.99]=NA
	counts=apply(octo,3,function(x) table(factor(x,levels=1:8)))
	
	# interaval estimation
	print("intervals for double check plot")
	delta=0.00001; nc=20
	interv1=seq(scft_conv-delta*nc,scft_conv+delta*nc,delta)
	interv0=seq(0,0.1,0.001)
	
	if (!do_parallel) {
		funct0=funct1=NULL
		t0=proc.time()
		funct0=sapply(interv0,function(x) bruteforce.fast(x))
		t1=proc.time(); print(t1-t0)
		funct1=sapply(interv1,function(x) bruteforce.fast(x))
		t2=proc.time(); print(t2-t1)
	}

	if (do_parallel) {
		cl <- makeCluster(cores)
		clusterExport(cl, varlist=list("scft_conv","heatflux","scbl_prime","w_prime","scft","np","lev","bruteforce.fast"))
		registerDoParallel(cl)
		#t0=proc.time()
		funct0=foreach(x=interv0,.combine=c) %dopar% bruteforce.fast(x)
		#t1=proc.time(); print(t1-t0)
		funct1=foreach(x=interv1,.combine=c) %dopar% bruteforce.fast(x)
		#t2=proc.time(); print(t2-t1)
		stopCluster(cl)
	}
         
        figname=file.path(FIGDIR,paste0(padnt,"_",expcode,"_convergence_intervals.pdf")); hh=10
        pdf(file=figname,height=hh,width=hh*2,onefile=TRUE,family="Helvetica",fonts="Helvetica")	
	par(mfrow=c(1,2),cex.main=3,cex.lab=3,cex.axis=3,mar=c(5,5,5,6),oma=c(1,1,1,5))
	plot(interv0,funct0,type="l")
	abline(v=scft_conv,col="blue")
	abline(h=maxflux,col="blue")

	plot(interv1,funct1,type="l")
	abline(v=scft_conv,col="blue")
	abline(h=maxflux,col="blue")
	dev.off()

	#Shape Factor
	print("shape factor...")
	circ_list=NULL
	octs=c(4,5,7)
	if (!do_parallel) {
		for (noct in octs) {
			print(noct)
			#t0=proc.time()
			circ=sapply(zz,function(z) {shape.factor.fast(lon,lat,octo[,,z],noct,min.area=50000)})
			#t1=proc.time()
        		#print(t1-t0)
			circ_list[[paste0("circ_",noct)]]<-circ
		}
	}
	if (do_parallel) {
		cl <- makeCluster(cores)
		clusterEvalQ(cl, library(SDMTools))
		clusterExport(cl, varlist=list("lon","lat","octo","zz","shape.factor.fast"))
		registerDoParallel(cl)
			for (noct in octs) {
				t0=proc.time()
				circ=foreach(z=zz,.combine=c) %dopar% shape.factor.fast(lon,lat,octo[,,z],noct,min.area=50000)
				#t1=proc.time(); print(t1-t0)
				circ_list[[paste0("circ_",noct)]]<-circ
			}
		stopCluster(cl)
	}

	#profile computes
	t0=proc.time()
	varprofiles=c("massflux","wtheta","wqtot","octo","qtot","tl","w","tke","scbl","scft","qrad3d","ql","buoyancy")

	#new code with sapply, faster by 30%
	profile_list=NULL
	for (var in varprofiles) {
		print(var)
		profile=profile.octants.new(var,octo,counts)
		profile=cbind(profile,profile.mean(var))
		profile_list[[var]]<-profile
	}

	save(file=save_profiles,profile_list,circ_list,lev,zi,zbmn,zcmn)
        
	}
	# ------------------------------- #
	# ------------ FIGURES ---------- # 
	# ------------------------------- #


	#useful names, palettes, values...
	TOP=1200
	LL=4
	PALOCT=c("lightpink","yellow","orange2","red3","palegreen","darkgreen","navy","royalblue")
	kk=1:length(counts[,1])
	kkk=c("Entrainment shell","Mixing Updraft","Ascending shell","Updraft","Entrainment","Mixing shell","Downdraft","Subsiding shell")
	heights=c(1,0.95,0.8,0.5,0.2,0.05)
	lettering=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)")
	cex.letter=3

	# ------------------------------- #
	# Figures of shapefactor
	print("shapefactor..")
	figname=file.path(FIGDIR,paste0(padnt,"_",expcode,"_shape_factor.pdf")); hh=15
	pdf(file=figname,height=hh,width=hh*0.5,onefile=TRUE,family="Helvetica",fonts="Helvetica")
	par(cex.main=3,cex.lab=2,cex.axis=2,mar=c(5,5,5,5))
	XX=c(0,0.25)
	plot(circ_list[["circ_7"]],lev[zz],col="navy",type="l",lwd=LL,xlab="Circularity",ylab="Height (m)",ylim=c(0,TOP),
	     main="Updraft/Downdraft Circularity",xlim=XX)+grid()
	points(circ_list[["circ_4"]],lev[zz],col="red",type="l",lwd=LL)
	points(circ_list[["circ_5"]],lev[zz],col="navy",type="l",lty=3,lwd=LL)
	abline(h=zi,col="darkgreen")
	legend(min(XX),TOP,c(kkk[c(4,7)],"Downdraft+Entrainment"),col=c("red","navy","navy"),lwd=LL,lty=c(1,1,3),cex=2)
	dev.off()


	# ------------------------------- #
	# Figures of octants flux profiles 
	print("Octants field profiles..")
	figname=file.path(FIGDIR,paste0(padnt,"_",expcode,"_field_profiles.pdf")); hh=40
	pdf(file=figname,height=hh,width=hh,onefile=TRUE,family="Helvetica",fonts="Helvetica")
	par(mfrow=c(2,4),cex.main=4.5,cex.lab=3.5,cex.axis=3.5,mar=c(5,5,5,5))

	varoctant=c("octo","qtot","tl","scft","scbl","w","tke","buoyancy")
	
	for (var in varoctant) {
		print(var); PPP=PALOCT
		if (var=="octo") {var0=counts; XX=c(0,50); XN="Percentage (%)"; MM="Cluster frequency"}
		if (var=="w") {var0=w; XX=c(-1.3,1.3); XN="Velocity (m/s)"; MM="Average vertical velocity"}
		if (var=="tl") {var0=tl; XX=c(289,290.5); XN="Temperature (K)"; MM="Liquid Potential Temperature"}
		if (var=="qtot") {var0=qtot; XX=c(0.006,0.0095); XN="Mixing ratio (g/kg)"; MM="Total water"}
		if (var=="tke") {var0=tke; XX=c(0,2); XN="TKE (m^2/s^2)"; MM="Turbulent Kinetic Energy"}
		if (var=="scft") {var0=scft; XX=c(0,1); XN=""; MM="Tracer2"}
		if (var=="scbl") {var0=scbl; XX=c(0,1.5); XN=""; MM="Scalar"}
		if (var=="buoyancy") {var0=buoyancy; XX=c(-0.015,0.015); XN="m/s^2" ; MM="Buoyancy"}
	
		plot(1,1,type="n",ylim=c(0,TOP),main=MM,ylab="Height (m)",lwd=LL,xlim=XX,xlab=XN)
		if (var!="octo") {points(profile_list[[var]][,9],lev,type="l",lwd=LL,lty=3)}

		abline(h=c(heights*zi),col="gray60",lwd=2)
		abline(h=c(zbmn,zcmn),col="darkgreen",lty=3,lwd=2)
		abline(h=zi,col="darkgreen",lwd=2)
		grid()
		mtext(lettering[which(var==varoctant)],line=1.5,adj=0.02,cex=cex.letter)
	
		for (noct in kk) {
			points(profile_list[[var]][,noct],lev,type="l",lwd=LL,col=PPP[noct])
		}
	abline(v=0)
	legend(min(XX),TOP,paste(kk,kkk),col=PALOCT[kk],cex=3,lwd=LL)
	}

	dev.off()

	# ------------------------------- #
	# Figures of octants flux profiles 
	print("Octants flux profiles..")
	figname=file.path(FIGDIR,paste0(padnt,"_",expcode,"_fluxes_profiles.pdf")); hh=18
	pdf(file=figname,height=hh,width=hh*5/3,onefile=TRUE,family="Helvetica",fonts="Helvetica")
	par(mfrow=c(1,3),cex.main=4.5,cex.lab=3.5,cex.axis=3.5,mar=c(5,5,5,5))

	varoctant=c("massflux","wtheta","wqtot")
	for (var in varoctant) {
		print(var); PPP=PALOCT
        	if (var=="wtheta") {var0=w_prime*tl_prime; XX=c(-0.035,0.035); XN="Heat flux (K*m/s)"; MM="Turbulent heat flux"}
        	if (var=="massflux") {var0=massflux; XX=c(-0.4,0.4); XN=""; MM="Turbulent mass flux"; XN="Mass flux (kg/(m^2*s))"}
        	if (var=="wqtot") {var0=w_prime*qtot_prime; XX=c(-0.6,0.6)*10^-4; XN="Moisture flux (g/kg*m/s)"; MM="Turbulent total water flux" }

        	plot(1,1,type="n",ylim=c(0,TOP),main=MM,ylab="Height (m)",lwd=LL,xlim=XX,xlab=XN)
        	points(profile_list[[var]][,9],lev,type="l",lwd=LL,lty=3)
	
        	abline(h=c(heights*zi),col="gray60",lwd=2)
        	abline(h=c(zbmn,zcmn),col="darkgreen",lty=3,lwd=2)
        	abline(h=zi,col="darkgreen",lwd=2)
        	grid()
		mtext(lettering[which(var==varoctant)],line=1.5,adj=0.02,cex=cex.letter)

		for (noct in kk) {
			points(profile_list[[var]][,noct],lev,type="l",lwd=LL,col=PPP[noct])
		}

		abline(v=0)
		legend(min(XX),TOP,paste(kk,kkk),col=PALOCT[kk],cex=3,lwd=LL)
		}

	dev.off()



if (do_sink) {sink()}
}



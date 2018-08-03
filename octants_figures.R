#program for octants and figures for dycoms2_rf01 stratocumulus convection
library(fields)
CODEDIR="/home/paolo/StratoCostaud"
source(file.path(CODEDIR,"les_routines.R"))

#setup flags
do_loading=T
do_computation=T
do_figures=T

#timestep 
nt=31
expcode="tecs00"


# ------------------------------- #
# ------- START PROGRAM --------- # 
# ------------------------------- #
model="uclales2"
exptype="dycoms2_rf01"
expname=paste(exptype,expcode,sep="_")
print(paste("Expname:",expname))

# configuration
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

# figure dirs
FIGDIR=file.path(FIGBASE,expname,"evolution",nt)
dir.create(FIGDIR,recursive=T)


# experiment-dependent properties
if (exptype=="dycoms2_rf01") {
	        u0=7; v0=-5.5
}

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
	
}

if (do_computation) {

	#PBL height defition
	zi=lev[which.max(diff(apply(tl,c(3),mean,na.rm=T)))] #end of the inversion layer, maximum gradient of T

	#primes
        print("eddies...")
        for (var in c("w","scbl","scft","tl","qtot")) {
                assign(paste0(var,"_prime"),eddy.mean(get(var)))
        }

	#select boundary layer points
	zz=1:whicher(lev,zi)

	# load from bruteforce the threshold to define scft_star	
	BRUTEDIR=file.path(FIGBASE,expname,"evolution",nt)
	savefile=file.path(BRUTEDIR,paste(padnt,expcode,"bruteforce_scft_conv",sep="_"))
	load(savefile)

	# computing octants
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

	# quadrant
	print("quadrant")
	quadra=w_prime*NA
	quadra[w_prime>0 & scbl_prime>0]=3
	quadra[w_prime<0 & scbl_prime>0]=7
	quadra[w_prime>0 & scbl_prime<0]=4
	quadra[w_prime<0 & scbl_prime<0]=8

	# cleans quadrants
	quadra[scft>0.99]=NA
	counts_quadra=apply(quadra,3,function(x) table(factor(x,levels=1:9)))

}
# ------------------------------- #
# ------------ FIGURES ---------- # 
# ------------------------------- #

if (do_figures) {

	# useful names, palettes, values...
	TOP=1200
	lpal=50
	PALREDBLUE=colorRampPalette(c("darkblue","blue","white","red","darkred"))(lpal)
	PALTIM=tim.colors(lpal)
	PALBLUEWHITE=colorRampPalette(c("darkblue","white"))(lpal)
	nocts=length(octname)
	kk=1:nocts
	heights=c(1,0.95,0.8,0.5,0.2,0.05)
	cex.letter=4

	# ------------------------------- #
	#Figures of horzontal sections
	#print("Horizontal sections for variables")
	#png(filename=paste0(FIGDIR,exp,"_var2d_sections.png"),height=4000,width=2550)
	#par(mfrow=c(8,5),cex.main=4,cex.lab=3,cex.axis=3,oma=c(1,1,1,3),mar=c(5,5,5,6))
	#
	#var2d=c("w","tracer","scft","scbl","octo")
	#heights2=c(1,0.95,0.9,0.8,0.5,0.2,0.1,0.05)
	#for (height in c(rev(heights2)))
	#        {
	#        vzi=whicher(lev,height*zi)
	#        for (var in var2d)
	#        {
	#        octo0=get(var)
	#        if (var!="octo")
	#        {image.plot(lon,lat,octo0[,,vzi],main=paste0(var," at ",round(height*zi,1),"m (",height,"*zi)"),xlab="X (m)",ylab="Y (m)",col=PALTIM)}
	#        else
	#        {image.plot(lon,lat,octo0[,,vzi],main=paste0("Octants at ",round(height*zi,1),"m (",height,"*zi)"),xlab="X (m)",ylab="Y (m)",col=PALOCT,zlim=c(1,nocts))}
	#        }
	#        }
	#
	#dev.off()

	# ------------------------------- #
	#Figures of vertical sections
	print("Vertical crossection...")
	varsection=c("ql","scbl","scft","w","quadra","octo")
	xx=whicher(lon,0)
	h0=1500

	png(filename=file.path(FIGDIR,paste0(expcode,"_cross_section.png")),height=h0*1.2,width=2*h0)
	par(mfrow=c(3,2),cex.main=5.5,cex.lab=3.5,cex.axis=4,mar=c(7,11,5,16),oma=c(1,0,7,1),mgp=c(6,2,0))
	for (var in varsection) {
		for (axis in c("lon")) {

			xx=whicher(get(axis),0); factor=1; mline=8
			if (axis=="lon") {var0=get(var)[xx,,]; XX="Y (m)"}
			if (axis=="lat") {var0=get(var)[,xx,]; XX="X (m)"}
			if (var=="octo") {PAL=PALOCT; zzz=c(1,8); mmm=""; MM="Octants"}
			if (var=="quadra") {PAL=PALOCT; zzz=c(1,8); mmm=""; MM="Quadrants"}
			if (var=="scbl") {zzz=c(0,1.4); PAL=rev(tim.colors()); MM=expression(S[BL]); mmm=""}
			if (var=="scft") {zzz=c(0,1); PAL=c(PALTIM,rep(PALTIM[lpal],lpal*9)); MM=expression(S[FT]); mmm=""}
			if (var=="w") {zzz=c(-3.6,3.6); PAL=PALREDBLUE; MM="Vertical velocity"; mmm="m/s"; mline=5}
			if (var=="ql") {zzz=c(0,0.5); PAL=PALBLUEWHITE; MM="Liquid water mixing ratio"; mmm=expression(paste(10^-3,"g/kg")); factor=1000 }
			
			#image(get(axis),lev,var0*factor,main=paste0(MM,": Vertical Cross-Section at lon=",lon[xx],"m"),ylab="",col=PAL,xlab="",zlim=zzz,ylim=c(0,TOP))
			image(get(axis),lev,var0*factor,main=bquote(bold(.(MM[[1]]))),ylab="",col=PAL,xlab="",zlim=zzz,ylim=c(0,TOP))
			#rect(-1700,200,-1000,900,lwd=8,border="gray80",lty=3)
			mtext("Altitude (m)",side=2,line=6,cex=3)
			mtext(XX,side=1,line=5,cex=3)
			if (var=="octo") {mtext(paste0("Vertical Cross-Section at X=",lon[xx],"m"),side=3,outer=T,line=2,cex=5,font=2)}
			#filled.contour3(get(axis),lev,var0,main=paste(MM,"Vertical Cross-Section at",lon[xx]),ylab="Altitude (m)",color.palette=PAL,xlab=XX,levels=seq(zzz[1],zzz[2],lpal),ylim=c(0,TOP))
			mtext(lettering[which(var==varsection)],line=1.5,adj=0.02,cex=cex.letter)
		}

		if (var=="octo" | var=="quadra") {
			image.plot(var0,col=PAL,zlim=zzz,legend.only=T,legend.width=5,
				   legend.args=list(side=4,line=6,cex=2.5,text=mmm),
				   axis.args=list(cex.axis=3.5,at=1:8),smallplot=c(0.92,0.94,0.15,0.9))
		}
		else {
			image.plot(var0,col=PAL,zlim=zzz,legend.only=T,legend.width=5,
				   legend.args=list(side=4,line=mline,cex=2.5,text=mmm),
				   axis.args=list(cex.axis=3.5),smallplot=c(0.92,0.94,0.15,0.9))
		}
	}
	dev.off()

	#------------------------------- #
	#Figures of horzontal sections
	print("Horizontal sections for variables")
	h0=1800
	png(filename=file.path(FIGDIR,paste0(expcode,"_horizontal_sections.png")),height=h0,width=32.6/20*h0)
	par(mfrow=c(2,3),cex.main=6,cex.lab=4,cex.axis=3.5,oma=c(1,1,7,1),mar=c(6,10,5,13),mgp=c(5,2,0))

	var2d=c("ql","w","quadra","scbl","scft","octo")
	heights2=c(0.9)
	for (height in c(rev(heights2))) {
        	vzi=whicher(lev,height*zi)
        	for (var in var2d) {
        		octo0=get(var); factor=1; mline=9
			if (var=="octo") {PAL=PALOCT; zzz=c(1,8); mmm=""; MM="Octants"}
        		if (var=="quadra") {PAL=PALOCT; zzz=c(1,8); mmm=""; MM="Quadrants"}
        		if (var=="scbl") {zzz=c(0.3,0.9); PAL=rev(tim.colors()); MM=expression(S[BL]); mmm=""}
        		if (var=="scft") {zzz=c(0,0.15); PAL=PALTIM; MM=expression(S[FT]); mmm=""}
        		if (var=="w") {zzz=c(-3.6,3.6); PAL=PALREDBLUE; MM="Vertical velocity"; mmm="m/s"; mline=5}
        		if (var=="ql") {zzz=c(0,0.40); PAL=PALBLUEWHITE; MM="Liquid water mixing ratio"; mmm=expression(paste(10^-3,"g/kg")); factor=1000 }


			image(lon,lat,factor*octo0[,,vzi],main=bquote(bold(.(MM[[1]]))),xlab="X (m)",ylab="Y (m)",col=PAL,zlim=zzz,asp=1)
			mtext(lettering[which(var==var2d)],line=1.5,adj=0.02,cex=cex.letter)
			if (var=="octo") {
				mtext(paste0("Horizontal Section at ",round(height*zi,1),"m (",height,"*zi)"),side=3,outer=T,line=2,cex=5,font=2)
			}
	
			if (var=="octo" | var=="quadra") {
				image.plot(var0,col=PAL,zlim=zzz,legend.only=T,legend.width=5,
					   legend.args=list(side=4,line=6,cex=2.5,text=mmm),
					   axis.args=list(cex.axis=4,at=1:8),smallplot=c(0.90,0.92,0.1,0.9))
			} else {
				image.plot(var0,col=PAL,zlim=zzz,legend.only=T,legend.width=5,
					   legend.args=list(side=4,line=mline,cex=2.5,text=mmm),
					   axis.args=list(cex.axis=4),smallplot=c(0.90,0.92,0.1,0.9))
			}
		}
        }

	dev.off()



	# ------------------------------- #
	#Figures of octants horizonal sections
	print("Horizontal sections...")

	h0=2000
	png(filename=file.path(FIGDIR,paste0(expcode,"_octants_sections.png")),height=h0,width=1.51*h0)
	par(mfrow=c(2,3),cex.main=6,cex.lab=4,cex.axis=4,oma=c(1,1,1,3),mar=c(8,8,5,5),mgp=c(5,1.5,0))

	for (height in c(rev(heights))) {
		vzi=whicher(lev,height*zi)
	
		image(lon,lat,octo[,,vzi],main=paste0("Octants at ",round(height*zi,1),"m (",height,"*zi)"),
		      xlab="X (m)",ylab="Y (m)",col=PALOCT,zlim=c(1,nocts),asp=1,useRaster=T)
		mtext(lettering[which(height==rev(heights))],line=1.5,adj=0.02,cex=cex.letter)
	}
	image.plot(octo[,,vzi],col=PALOCT,zlim=c(1,nocts),legend.only=T,legend.width=5,axis.args=list(cex.axis=3.5,at=1:8))

	dev.off()

	# ------------------------------- #
	# Figures of octants scatterplot 
	print("octants scatterplot...")

	#rx=sample(1:length(lon),length(lon)/2); ry=sample(1:length(lat),length(lat)/2)
	rx=1:(length(lon)/2); ry=1:(length(lat)/2)
	CCC=0.1 
	var1="tl"
	var2="qtot"
	bdcol=c("black","black","black","white","black","white","white","white")
	png(filename=file.path(FIGDIR,paste0(expcode,"_",var1,"_",var2,"_scatterplot2.png")),height=2000,width=3100)
	#figname=paste0(FIGDIR,exp,"_",var1,"_",var2,"_scatterplot.pdf"); hh=40
	#pdf(file=figname,height=hh,width=hh*31/20,onefile=TRUE,family="Helvetica",fonts="Helvetica")
	par(mfrow=c(2,3),cex.main=6,cex.lab=4,cex.axis=4,mar=c(8,10,7,3),mgp=c(6,2,0),pty="s")

	for (height in c(rev(heights))) {
		XXX=c(289,290); MM1="Liquid Potential Temperature (K)"; 
		YYY=c(0.0084,0.0096); MM2="Total water mixing ratio (kg/kg)";
		vzi=whicher(lev,height*zi); hzi=round(height*zi,1)
		if (height==0.95) {XXX=c(289,292); YYY=c(0.0075,0.0096)}
		if (height==1) {XXX=c(289,300); YYY=c(0.0015,0.0096)}
		tt=substitute(paste(theta[l],"-",q[tot]," scatterplot at ",hz,"m (",h,"*zi)"),list(h=height,hz=hzi))

		#plot(get(var1)[rx,ry,vzi],get(var2)[rx,ry,vzi],col=PALOCT[octo[rx,ry,vzi]],main=tt,xlab="",ylab="",xlim=XXX,ylim=YYY,cex=CCC)
		plot(get(var1)[,,vzi],get(var2)[,,vzi],col=PALOCT[octo[,,vzi]],main=tt,xlab="",ylab="",xlim=XXX,ylim=YYY,cex=CCC)
		grid()
		#plot(get(var1)[,,vzi],get(var2)[,,vzi],col=PALOCT[octo[,,vzi]],main=tt,xlab=MM1,ylab=MM2,xlim=XXX,ylim=YYY,cex=CCC)
		mtext(MM2,side=2,line=6,cex=3)
        	mtext(MM1,side=1,line=6,cex=3)
		abline(h=mean(get(var2)[,,vzi]))
		abline(v=mean(get(var1)[,,vzi]))
		legend(max(XXX)-(max(XXX)-min(XXX))/2.5,max(YYY),paste(kk,octname),col=PALOCT[kk],cex=4,pch=19,bg="white")
		#points(c(tl_ml,tl_ft),c(qtot_ml,qtot_ft),type="l",cex=2,lty=3)
		mtext(lettering[which(height==rev(heights))],line=1.5,adj=0,cex=cex.letter)
		
		for (nn in 1:nocts) {
			v1=get(var1)[,,vzi]; v2=get(var2)[,,vzi]
			points(mean(v1[octo[,,vzi]==nn],na.rm=T),mean(v2[octo[,,vzi]==nn],na.rm=T),bg=PALOCT[nn],pch=23,cex=10,col=bdcol[nn])
			#v1=get(var1); v2=get(var2)
			#points(mean(v1[octo[,,vzi]==nn],na.rm=T),mean(v2[octo[,,vzi]==nn],na.rm=T),col="black",bg=PALOCT[nn],pch=23,cex=12)
		}
	}
	dev.off()


}

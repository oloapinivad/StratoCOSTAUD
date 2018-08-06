#averaging profiles
#rm(list=ls())
CODEDIR="/home/paolo/StratoCostaud"
source(file.path(CODEDIR,"config.R"))
source(file.path(CODEDIR,"les_routines.R"))

# CONFIGURATION 
expcode="tecs00"
exptype="dycoms2_rf01"
expname=paste(exptype,expcode,sep="_")
ntimes=1:31

print(paste("Expname:",expname))

FIGDIR=file.path(FIGBASE,expname)
dir.create(FIGDIR,recursive=T)

#preloading
load(file.path(FIGBASE,expname,"evolution","01",paste("01",expcode,"profiles",sep="_")))
varprofile=names(profile_list)
varcirc=names(circ_list)
vartime=c("zi","scft_conv","scft_zi")
for (var in vartime)	{assign(paste0("full_",var),1:length(ntimes)*NA)}
for (var in varprofile) {assign(paste0("full_profile_",var),array(NA,dim=c(length(lev),9,length(ntimes))))}
for (var in varcirc)	{assign(paste0("full_",var),array(NA,dim=c(length(lev),length(ntimes))))}

for (nt in ntimes) {

	# load profiles
	padnt=sprintf("%02d",nt)
	LOADDIR=file.path(FIGBASE,expname,"evolution",padnt)
	save_file=file.path(LOADDIR,paste(padnt,expcode,"bruteforce_scft_conv",sep="_"))
	save_profiles=file.path(LOADDIR,paste(padnt,expcode,"profiles",sep="_"))
	load(save_file)
	load(save_profiles)

	pnt=nt-min(ntimes)+1

	for (var in vartime) {
		full=get(paste0("full_",var))
		full[pnt]=get(var)
		assign(paste0("full_",var),full)
	}

	for (var in varprofile) {
		full_profile=get(paste0("full_profile_",var))
		full_profile[,,pnt]=profile_list[[var]]
		assign(paste0("full_profile_",var),full_profile)
	}

	for (var in varcirc) {
		full_circ=get(paste0("full_",var))
		ll=length(circ_list[[var]])
		full_circ[1:ll,pnt]=circ_list[[var]]
		assign(paste0("full_",var),full_circ)
	}
}

#compute mean
for (var in varprofile) {
	assign(paste0("mean_profile_",var),rowMeans(get(paste0("full_profile_",var)),dim=2,na.rm=T))
}

for (var in varcirc) {
	assign(paste0("mean_",var),apply(get(paste0("full_",var)),1,mean,na.rm=T))
}

#useful names, palettes, values...
TOP=1200
kk=1:length(kkk)
heights=c(1,0.95,0.8,0.5,0.2,0.05)
cex.letter=3
LL=6
zi_mean=mean(zi)

# ------------------------------- #
# Figures of time evolution of entrainment
#print("timevolution..")
#LL=4
#png(filename=paste0(FIGDIR,"mean_",exp,"_fluxes_evolution.png"),height=2000,width=2000)
#par(cex.main=4,cex.lab=3,cex.axis=3,mar=c(5,5,5,5),mfrow=c(5,1))
#plot(ntimes*120,full_scft_conv,type="l",lwd=LL,main="Threshold (solid) and Zi-average tracer (dashed) and 800-850m averaged tracer (green)",xlab="time (sec)",ylim=c(0.01,0.03),ylab="Tracer concentration [0-1]")
#points(ntimes*120,full_scft_zi,type="l",lwd=LL,lty=2,col="red")
#points(ntimes*120,full_tr2ml,type="l",lwd=LL,lty=2,col="darkgreen")
#text(120,0.03,paste("Cor.= ",round(cor(full_scft_zi,full_scft_conv),2)),cex=3,pos=4)
#grid()
#plot(ntimes*120,full_zi,type="l",lwd=LL,main="Zi Evolution",xlab="time (sec)",ylim=c(860,880),col="red",ylab="Height (m)")
#grid()
#plot(ntimes*120,cloud_flux,type="l",lty=2,lwd=LL,main="Cloud base integral heat flux (600-650m)",xlab="time (sec)",ylim=c(-0.3,0.3),ylab="K*m^2/s")
#text(120,0.03,paste("Cor (with threshold)= ",round(cor(cloud_flux,full_scft_conv),2)),cex=3,pos=4)
#grid()
#plot(ntimes*120,abs(entrain_flux),type="l",col=PALOCT[5],lwd=LL,ylim=c(2,4),main="Integral of Entrainment and Downdraft flux",xlab="time (sec)",ylab="K*m^2/s")
#points(ntimes*120,abs(down_flux),type="l",col=PALOCT[7],lwd=LL)
##names=c("20x20x1A (256)","20x20x5 (256)","10x10x5 (256)", "20x20x5 (256) - FullRad","20x20x5 (256) - RightRad","20x20x5 (256) - SST")
#text(120,2.5,paste("Cor.= ",round(cor(down_flux,entrain_flux),2)),cex=3,pos=4)
#text(120,2.3,paste("Cor. w. cb flux= ",round(cor(down_flux,cloud_flux),2)),cex=3,pos=4,col=PALOCT[7])
##text(120,2.1,paste("Cor. w. cb flux= ",round(cor(entrain_flux,cloud_flux),2)),cex=3,pos=4,col=PALOCT[5])
#ccf(full_tr2ml,cloud_flux,main="Lagged correlation 800-850m average tracer (leading) vs. Cloud heat flux")
#grid()
#dev.off()

# ------------------------------- #
# Figures of shapefactor
print("shapefactor..")
figname=file.path(FIGDIR,paste0("mean_",expcode,"_shape_factor.pdf")); hh=15
pdf(file=figname,height=hh,width=hh*0.5,onefile=TRUE,family="Helvetica",fonts="Helvetica")
par(cex.main=3,cex.lab=2,cex.axis=2,mar=c(5,5,5,5))
XX=c(0,0.25)
plot(mean_circ_7,lev,col="navy",type="l",lwd=LL,xlab="Circularity",ylab="Height (m)",ylim=c(0,TOP),main="Updraft/Downdraft Circularity",xlim=XX)+grid()
points(mean_circ_4,lev,col="red",type="l",lwd=LL)
points(mean_circ_5,lev,col="navy",type="l",lty=3,lwd=LL)
abline(h=zi,col="darkgreen")
legend(min(XX),TOP,c(kkk[c(4,7)],"Downdraft+Entrainment"),col=c("red","navy","navy"),lwd=LL,lty=c(1,1,3),cex=2)
dev.off()

#weighted.mean(mean_circ_7[1:length(diff(lev))],diff(lev),na.rm=T)

# ------------------------------- #
# Figures of octants flux profiles 
print("Octants field profiles..")
#png(filename=paste0(FIGDIR,"mean_",exp,"_field_profiles.png"),height=3000,width=3500)
figname=file.path(FIGDIR,paste0("mean_",expcode,"_field_profiles.pdf")); hh=40
pdf(file=figname,height=hh,width=hh,onefile=TRUE,family="Helvetica",fonts="Helvetica")
par(mfrow=c(2,4),cex.main=4.5,cex.lab=4,cex.axis=4,mar=c(9,11,5,3),mgp=c(7,2,0))

varoctant=c("octo","qtot","tl","ql","w","buoyancy","tke","qrad3d")
for (var in varoctant) {

        print(var); PPP=PALOCT; factor=1
	profile=get(paste0("mean_profile_",var))
        if (var=="octo") {XX=c(0,50); XN="Percentage (%)"; MM="Cluster frequency"}
        if (var=="w") {XX=c(-1.3,1.3); XN="Velocity (m/s)"; MM="Average vertical velocity"}
        if (var=="qrad3d") {XX=c(-10,2); XN="Cooling rate (K/hour)"; MM="Radiative cooling" }
        if (var=="tl") {XX=c(289,290.5); XN="Temperature (K)"; MM="Liquid Potential Temperature"}
        if (var=="qtot") {XX=c(0.006,0.0095); XN="Mixing ratio (kg/kg)"; MM="Total water"}
        if (var=="tke") {XX=c(0,2); XN=expression(paste("TKE (",m^2/s^2,")")); MM="Turbulent Kinetic Energy"}
        if (var=="scft") {XX=c(0,1); XN=""; MM="Tracer2"}
        if (var=="scbl") {XX=c(0,1.5); XN=""; MM="Scalar"}
        if (var=="buoyancy") {XX=c(-0.020,0.020); XN=expression(m/s^2) ; MM="Buoyancy"}
	if (var=="ql") {XX=c(0,0.6); XN=expression(paste("Mixing ratio (",10^-3,g/kg,")")) ; MM="Liquid water"; factor=10^3}

        plot(1,1,type="n",ylim=c(0,TOP),main=MM,ylab="Height (m)",lwd=LL,xlim=XX,xlab=XN)
	#mtext("Height (m)",side=2,line=5,cex=3)
        #mtext(XN,side=1,line=5,cex=3)
        points(profile[,9]*factor,lev,type="l",lwd=LL,lty=3)

        #abline(h=c(heights*zi),col="gray60",lwd=2)
        #abline(h=c(zbmn,zcmn),col="darkgreen",lty=3,lwd=2)
        abline(h=zi,col="darkgreen",lwd=2)
        grid()
        mtext(lettering[which(var==varoctant)],line=1.5,adj=0,cex=cex.letter)

        for (noct in kk) {
                points(profile[,noct]*factor,lev,type="l",lwd=LL,col=PPP[noct])
        }
        abline(v=0)
        legend(min(XX),TOP,paste(kk,kkk),col=PALOCT[kk],cex=3,lwd=LL,bg="white")
}
#
#XX=c(0,0.25)
#plot(mean_circ_7,lev,col="navy",type="l",lwd=LL,xlab="Circularity",ylab="Height (m)",ylim=c(0,TOP),main="Updraft/Downdraft Circularity",xlim=XX)+grid()
#points(mean_circ_4,lev,col="red",type="l",lwd=LL)
#points(mean_circ_5,lev,col="navy",type="l",lty=3,lwd=LL)
#abline(h=zi,col="darkgreen")
#legend(min(XX),TOP,c(kkk[c(4,7)],"Downdraft+Entrainment"),col=c("red","navy","navy"),lwd=LL,lty=c(1,1,3),cex=2,bg="white")

dev.off()

# ------------------------------- #
# Figures of octants flux profiles 
print("Octants flux profiles..")
figname=file.path(FIGDIR,paste0("mean_",expcode,"_fluxes_profiles.pdf")); hh=18
pdf(file=figname,height=hh,width=hh*5/3,onefile=TRUE,family="Helvetica",fonts="Helvetica")
#png(filename=paste0(FIGDIR,"mean_",exp,"_fluxes_profiles.png"),height=1500,width=2500)
par(mfrow=c(1,3),cex.main=5,cex.lab=3.5,cex.axis=3.5,mar=c(10,10,5,5),mgp=c(5,2,0))

varoctant=c("massflux","wtheta","wqtot")
for (var in varoctant)
{
        print(var); PPP=PALOCT
	profile=get(paste0("mean_profile_",var))
        if (var=="wtheta") {XX=c(-0.035,0.035); XN=expression(paste("Heat flux (",K*m*s^-1,")")); MM="Turbulent heat flux"}
        if (var=="massflux") {XX=c(-0.4,0.4); XN=""; MM="Turbulent mass flux"; XN=expression(paste("Mass flux (",kg*m^-2*s^-1,")"))}
        if (var=="wqtot") {XX=c(-0.6,0.6)*10^-4; XN=expression(paste("Moisture flux (",kg*kg^-1*m*s^-1,")")); MM="Turbulent total water flux" }

        plot(1,1,type="n",ylim=c(0,TOP),main=MM,ylab="",lwd=LL,xlim=XX,xlab="")
        points(profile[,9],lev,type="l",lwd=LL,lty=3)
	mtext("Height (m)",side=2,line=5.5,cex=2.8)
        mtext(XN,side=1,line=7,cex=2.8)

        #abline(h=c(heights*zi),col="gray60",lwd=2)
        #abline(h=c(zbmn,zcmn),col="darkgreen",lty=3,lwd=2)
        abline(h=zi,col="darkgreen",lwd=2)
        grid()
        mtext(lettering[which(var==varoctant)],line=1.5,adj=0.02,cex=cex.letter)

        for (noct in kk)
                {
                points(profile[,noct],lev,type="l",lwd=LL,col=PPP[noct])
                }

        abline(v=0)
        legend(min(XX),TOP,paste(kk,kkk),col=PALOCT[kk],cex=3,lwd=LL,bg="white")
        }

dev.off()


#print some average quantity
#for (n in c(4,5,7)) {
#	avg=weighted.mean(get(paste0("mean_circ_",n))[1:length(diff(lev))],diff(lev),na.rm=T)
#	print(paste("Mean Circularity for octant",n,"is",round(avg,5)))
#}

#intflux=array(NA,dim=c(9,3))
#for (var in varoctant) {
#	print(var)
#	profile=get(paste0("mean_profile_",var))
#	for (noct in c(kk,9)) {
#		totflux=sum(profile[1:length(diff(lev)),noct]*diff(lev),na.rm=T)
#		intflux[noct,which(var==varoctant)]=totflux
#		print(paste("Integrated",var,"for octant",noct,"is:",round(totflux,3)))
#	}
#}

#print("Octants momentum profiles..")
#png(filename=paste0(FIGDIR,"mean_",exp,"_momentum_profiles.png"),height=1500,width=3000)
#par(mfrow=c(1,6),cex.main=3,cex.lab=2.5,cex.axis=2.5,mar=c(5,5,5,5))
#
#LL=4
#varoctant=c("buoyancy","advection","pressure","tendency","diffusion","residual")
#nocts=kk
#PPP=PALOCT
#klegend=paste(kk,kkk)
#XXXX=c(-0.005,0.005)
#YYYY=c(0,TOP)
#
#for (var in varoctant)
#{
#        print(var)
#	profile=get(paste0("mean_profile_",var))
#        XX=c(-0.005,0.005); XN="m/s^2"; MM=var
#
#        plot(1,1,type="n",ylim=YYYY,main=MM,ylab="Height (m)",lwd=LL,xlim=XX,xlab=XN)
##        points(profile[,9],lev,type="l",lwd=LL,lty=3)
#
        #abline(h=c(zbmn,zcmn),col="darkgreen",lty=3,lwd=2)
#        abline(h=zi,col="darkgreen",lwd=2)
#        grid()
#        mtext(lettering[which(var==varoctant)],line=1.5,adj=0.02,cex=cex.letter)
#
#for (noct in kk)
##	{
#	points(profile[,noct],lev,type="l",lwd=LL,col=PPP[noct])
##	}
#
#abline(v=0)
##legend(min(XX),TOP,klegend,col=PPP[nocts],cex=2,lwd=LL)
#}
#

#dev.off()
#



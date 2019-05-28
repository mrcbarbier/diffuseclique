#process data downloaded from Cedar Creek Natural Reserve Database
data=read.delim("e120_Plant aboveground biomass data.txt",header=T)

dim(data)
#remove data==NA
data=data[!is.na(data$Biomass..g.m2.),]
dim(data)
#remove data==0
data=data[!(data$Biomass..g.m2.==0),]

#remove duplicated data entries
potdoublon=c()
for(i in 1:(dim(data)[1]-1))
if((data[i,"Biomass..g.m2."]==data[i+1,"Biomass..g.m2."])&
	(data[i,"Species"]==data[i+1,"Species"]))
potdoublon=c(potdoublon,i)

for(k in potdoublon) print(data[k+(0:1),])

#remove doublons
data=data[-c(18112,18267),]

#correct errors in dates
data[(data$Plot==224)&(data$Year==2006),"Month"]=7
# for other errors, pool july and august dates=78 & differentiate june=6
data[(data$Month==7),"Month"]=78
data[(data$Month==8),"Month"]=78

names(data)
attach(data)
table(data$Species)

# reformat the data: one column per identified species
ReFormattedData=c()
ListAll=names(table(Species))
Bzzz=rep(0, length(ListAll))
names(Bzzz)=ListAll
length(ListAll)

for(yr in unique(Year))
for(mt in unique(Month))
for(pl in unique(Plot))
for(st in unique(Strip)[-5])
for(sst in unique(Substrip))
{
gz=data[(!is.na(Strip))&(Year==yr)&(Month==mt)&(Plot==pl)&(Strip==st)&(Substrip==sst),]
if(dim(gz)[1]>0)
	{
	Newline=c(gz[1,-c(8,36,37)],Bzzz)
	for (k in 1:(dim(gz)[1]))
	Newline[as.character(gz[k,"Species"])]=as.numeric(Newline[as.character(gz[k,"Species"])])+gz[k,"Biomass..g.m2."]
	ReFormattedData=rbind(ReFormattedData,Newline)
	}
}

dim(ReFormattedData)
Tot=rowSums(matrix(as.numeric(ReFormattedData[,35:210]), nrow=dim(ReFormattedData)[1],
	ncol=dim(ReFormattedData[,35:210])[2]))

ReFormattedData=cbind(ReFormattedData, Tot)
ReFormattedData=as.matrix(ReFormattedData)
dim(ReFormattedData)
dimnames(ReFormattedData)

detach(data)
########################################################
attach(ReFormattedData)

# remove june measurements
table(Month)
Sdata=ReFormattedData[Month!=6,]

# remove unsorted samples
Sdata1=Sdata[round(Sdata$Miscellaneous.litter+Sdata$Unsorted.biomass,3)!=round(Sdata$Tot,3),]
Sdata2=Sdata1[round(Sdata1$Miscellaneous.litter+Sdata1$Unsorted.Biomass,3)!=round(Sdata1$Tot,3),]
SortedData=Sdata2[round(Sdata2$Miscellaneous.litter+Sdata2$Green.matter..alive.,3)!=round(Sdata2$Tot,3),]

SortedData[((SortedData$NumSp==1)&(SortedData$Achmi==1)),c("Year","Plot","Achillea.millefolium.lanulosa.","Miscellaneous.litter","Green.matter..alive.","Tot")]
SortedData[c("4777","4778"),]
table(SortedData$Year,SortedData$Substrip)

#data over substrips should be averaged = keep both

subP=unique(SortedData[SortedData$Substrip>1,"Plot"])
nonsubP=unique(SortedData[SortedData$Substrip==1,"Plot"])
mean(SortedData[(is.element(SortedData$Plot,subP)&(SortedData$Year==2007)),"Tot"])
mean(SortedData[(is.element(SortedData$Plot,subP)&(SortedData$Year==2006)),"Tot"])
mean(SortedData[(is.element(SortedData$Plot,subP)&(SortedData$Year==2008)),"Tot"])
mean(SortedData[(is.element(SortedData$Plot,nonsubP)&(SortedData$Year==2007)),"Tot"])
mean(SortedData[(is.element(SortedData$Plot,nonsubP)&(SortedData$Year==2006)),"Tot"])
mean(SortedData[(is.element(SortedData$Plot,nonsubP)&(SortedData$Year==2008)),"Tot"])

sd(SortedData[(is.element(SortedData$Plot,subP)&(SortedData$Year==2007)),"Tot"])
sd(SortedData[(is.element(SortedData$Plot,subP)&(SortedData$Year==2006)),"Tot"])
sd(SortedData[(is.element(SortedData$Plot,subP)&(SortedData$Year==2008)),"Tot"])
sd(SortedData[(is.element(SortedData$Plot,nonsubP)&(SortedData$Year==2007)),"Tot"])
sd(SortedData[(is.element(SortedData$Plot,nonsubP)&(SortedData$Year==2006)),"Tot"])
sd(SortedData[(is.element(SortedData$Plot,nonsubP)&(SortedData$Year==2008)),"Tot"])


## combine Amorpha and Petalostemum into one compound species
## combine Monarda and Solidago into one compound species

AmoPet=(SortedData$Amocan|SortedData$Petpu)*1
table(AmoPet)

SimSorData1=SortedData[,c(2,4,8,17,18,20:28,30:34)]
SimSorData1=cbind(SimSorData1,AmoPet)
names(SimSorData1)[13]="MonSol"
SimSorData1=cbind(SimSorData1,SortedData[,c(35,37,44,56,83,108,116,117,120)])
Monarda.p.Solidago=SortedData[,134]+SortedData[,188]
SimSorData1=cbind(SimSorData1,Monarda.p.Solidago)
SimSorData1=cbind(SimSorData1,SortedData[,c(147,160,169,170,176,190)])
Amorpha.p.Petalostemum=SortedData[,43]+SortedData[,152]+SortedData[,153]+SortedData[,154]+SortedData[,155]
SimSorData1=cbind(SimSorData1,Amorpha.p.Petalostemum)
OtherSp=rowSums(SortedData[,c(36,38:42,45:55,57:82,84:107,109:115,118:119,
	121:133,135:146,148:151,156:159,161:168,171:175,177:187,189,191:210)])
SimSorData1=cbind(SimSorData1,OtherSp)
names(SimSorData1)
Chk=round(rowSums(SimSorData1[,21:38])-SortedData$Tot,2)
which(Chk>0)
dim(SimSorData1)
                
SimSorData1[((SimSorData1$NumSp==1)&(SimSorData1$Achmi==1)),c("Year","Plot","Achillea.millefolium.lanulosa.")]

names(SimSorData1)
# identify problematic species: trees (Quercus), Elymus & Agropyron that barely grew in monoculture
# plots where they took over need to be excluded from the analysis
sum(SimSorData1$Quercus.macrocarpa>50)
sum(SimSorData1$Quercus.ellipsoidalis>50)
sum(SimSorData1$Elymus.canadensis>50)
sum(SimSorData1$Agropyron.smithii>50)

# identify issue plots: thoses where these species took over
PbPlots=c()
PbPlots=unique(c(PbPlots, SimSorData1[(SimSorData1$Queel==1)&(SimSorData1$Quercus.ellipsoidalis>50),"Plot"]))
PbPlots=unique(c(PbPlots, SimSorData1[(SimSorData1$Quema==1)&(SimSorData1$Quercus.macrocarpa>50),"Plot"]))
PbPlots=unique(c(PbPlots, SimSorData1[(SimSorData1$Elyca==1)&(SimSorData1$Elymus.canadensis>50),"Plot"]))

SimSorData2=SimSorData1[!is.element(SimSorData1$Plot,PbPlots),]
write.csv(SimSorData2, file="BigBioNoPbPlots.csv",row.names = F)
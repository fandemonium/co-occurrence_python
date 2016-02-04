# This is a demo for running the co-occurrence analysis much, much faster

#make sure you have these libraries
library(Hmisc)
library(plyr)
library(reshape2)
library(igraph)
library(fdrtool)
# this is the demo data, take a look at it
dataset<-read.table("~/Desktop/pit_foaming_otus.txt",header=T,sep="\t",check.names=F)

# we are going to create a network per treatment
head(dataset[,1:10])

treatments<-as.vector(unique(dataset$Foaming_Status))
final_results<-data.frame()
for(i in 1:length(treatments)){
	#subset the data for a particular treatment
	temp<-subset(dataset, Foaming_Status==treatments[i])
	# making an object that has all the results in it (both rho and P values)
	results<-rcorr(as.matrix(temp[,-c(1:2)]),type="spearman")
	#make two seperate objects for p-value and correlation coefficients
	rhos<-results$r
	ps<-results$P
	# going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
	ps_melt<-na.omit(melt(ps))
	#creating a qvalue based on FDR
	ps_melt$qval<-fdrtool(ps_melt$value, statistic="pvalue", plot=F,verbose=F)$qval
	#making column names more relevant
	
	names(ps_melt)[3]<-"pval"
	# if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
	ps_sub<-subset(ps_melt, qval < 0.05)
	# now melting the rhos, note the similarity between ps_melt and rhos_melt
	rhos_melt<-na.omit(melt(rhos))
	names(rhos_melt)[3]<-"rho"
	#merging together and remove negative rhos
	merged<-merge(ps_sub,subset(rhos_melt, rho > 0),by=c("Var1","Var2"))
	merged$trt<-treatments[i]
	final_results<-rbind(final_results, merged)
	print(paste("finished ",treatments[i],sep=""))
}

# to make a network from the results
# we need to pass the relevant relationships (pairs across columns Var1 and Var2) to graph.edgelist as a matrix, note I am subsetting for a particular treatment within the arguments
head(final_results)
g<-graph.edgelist(as.matrix(subset(final_results, trt=="foaming")[,1:2]),directed=F)

#are relationships different
ggplot(final_results)+geom_density(aes(rho,fill=trt),alpha=0.5)+theme_bw(base_size=17)+theme(aspect.ratio=1)+scale_fill_manual(name="Treatment",values=c("red","black"))

# now we can calculate stats for the network
final_stats<-data.frame()
for(i in 1:length(unique(final_results$trt))){
	
	temp<-subset(final_results, trt==as.vector(unique(final_results$trt))[i])
	temp.graph<-(graph.edgelist(as.matrix(temp[,c(1,2)]),directed=FALSE))
	E(temp.graph)$weight<-temp$rho
	temp.graph<-simplify(temp.graph)
	stats<-data.frame(row.names((as.matrix(igraph::degree(temp.graph,normalized=TRUE)))),(as.matrix(igraph::degree(temp.graph,normalized=TRUE))),(as.matrix(igraph::betweenness(temp.graph))))
	names(stats)<-c("otus","norm_degree","betweenness")
	stats$trt<-as.vector(unique(final_results$trt))[i]
	stats$clustering_coeff<-igraph::transitivity(temp.graph,type="global")
	stats$clustering_coeff_rand<-igraph::transitivity(igraph::erdos.renyi.game(length(V(temp.graph)),length(E(temp.graph)),type="gnm"))
	stats$cluster_ratio<-stats$clustering_coeff/stats$clustering_coeff_rand
	final_stats<-rbind(final_stats,stats)
	print(paste("finished ",as.vector(unique(final_results$trt))[i],sep=""))
}

# is the structure changing?
ddply(final_stats,.(trt),summarise,mean(cluster_ratio))

head(arrange(subset(final_stats,trt=="foaming"), -norm_degree),10)
OTU 166 d:Bacteria,p:Firmicutes,c:Negativicutes,o:Selenomonadales,f:Acidaminococcaceae
OTU 2080 d:Bacteria,p:"Bacteroidetes"
OTU 1864 d:Bacteria

head(arrange(subset(final_stats,trt=="non-foaming"), -betweenness),10)
OTU 1864 d:Bacteria
OTU 1601  d:Bacteria,p:"Proteobacteria",c:Betaproteobacteria

# whats the relationship between these statistics?
ggplot(final_stats)+geom_point(aes(x=norm_degree,y=betweenness,color=trt),alpha=0.5)+scale_y_log10()+theme_bw(base_size=17)+labs(x="Normalized Degree",y="Betweenness")+theme(aspect.ratio=1)+scale_colour_manual(name="Treatment",values=c("red","black"))




# let's make plots with the most significant relationships
head(final_results)
dim(final_results)
hist(final_results$rho)
strong_results<-subset(final_results, rho >= 0.7)
dim(strong_results)
hist(strong_results$rho)

temp.graph<-(graph.edgelist(as.matrix(subset(strong_results, trt=="foaming")[,c(1,2)]),directed=FALSE))
E(temp.graph)$weight<-subset(strong_results, trt=="foaming")$rho
temp.graph<-simplify(temp.graph)
temp.graph

library(GGally)
library(intergraph)
gnet<-asNetwork(temp.graph)
df<-asDF(gnet)
head(df$vertexes)

ggnet(gnet, size=0, method="kamadakawaii")+geom_point()

# lets color by phylum
vs<-df$vertexes
phyla<-read.table("~/Desktop/otu_phylum.txt",sep="\t",check.names=F,header=T)
head(phyla)
head(vs)
vs_phyla<-merge(vs, phyla, by.x="vertex.names",by.y="otus")
vs_phyla<-arrange(vs_phyla,intergraph_id)
head(vs_phyla)
ggnet(gnet, size=0, method="kamadakawaii")+geom_point(aes(colour=vs_phyla$phylum))+scale_colour_manual(values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","grey50"))
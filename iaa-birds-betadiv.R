#Code for
#Prasetya et al.
#Birds and Barriers: present and past seas are dominant 
#correlates of avian turnover in the Indo-Australia Archipelago.

#load libraries
library(betapart)
library(readr)
library(ape)
library(phytools)
library(viridis)
library(ggplot2)
library(reshape)
library(vegan)
library(reshape2)
library(igraph)
library(dendextend)
library(wesanderson)
library(gridExtra)

# iaabird.sp <- read.csv("jetzspcommatrix.csv")
# iaabird.gen <- read.csv("jetzgencommatrix.csv")
# colnames(iaabird.gen)[1] <- "Abroscopus" #fix name
# iaabird.f <- read_delim("jetzfamcommatrix.csv", "\t", escape_double = FALSE, trim_ws = TRUE)

iaabird.sp <- read.csv("spcommmatrix.csv")
iaabird.sp <- iaabird.sp[,-1]
iaabird.gen <- read.csv("gencommatrix.csv")
iaabird.gen <- iaabird.gen[,-1]
iaabird.f <- read.csv("famcommatrix.csv")
iaabird.f <- iaabird.f[,-1]

#write rownames for iaabird
#"AUSTRALIA","MELANESIA","MALUKU","SULAWESI","PHILLIPPINES","LESSER SUNDA","BORNEO","JAVA","SUMATRA","SOUTHEAST ASIA"
rownames(iaabird.sp) <- c("AUS","NWG", "EMN", "MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA")
rownames(iaabird.gen) <- c("AUS","NWG", "EMN", "MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA")
rownames(iaabird.f) <- c("AUS","NWG", "EMN", "MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA")

#TREE CLEANING + CREATe
####TREES FOR SPECIES
#trees=read.nexus("output.nex")
#tree <- trees[[1]]

# splist <- read.csv("masterjetzsplist.csv")
# colnames(splist) <- "PHYLO_NAME"
# splist$PHYLO_NAME <- gsub(" ", "_", splist$PHYLO_NAME)
# #tree pruned to only selected species
# tips <- data.frame(name=tree$tip.label) #tips in original tree
# tips1 <- subset(tips, !tips$name%in%splist$PHYLO_NAME) #tips that we don't want 
# #RECHECK We are missing 326
# drop <- as.character(droplevels(tips1$name))
# 
# jetz_tree <- drop.tip(tree, drop) #tree with only what we want
# unique(jetz_tree$tip.label) #check same as tips we want to keep

#check <- subset(splist, !splist$PHYLO_NAME%in%jetz_tree$tip.label) #2 tips remain, NEW and extinct
#write.csv(check, "maisiespcheckmissingfromjetz.csv")

#GENUS TREE
#genus only tree
# tree <- read.newick("jetztreecleaned.newick")
# tips <- tree$tip.label
# genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
# 
# #drop all species but one of each
# ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
# 
# tree1<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
# tree1$tip.label
# 
# #Anser gone, adding Anser manually
# tree2 <- bind.tip(tree1, "Anser_indicus", where=which(tree$tip.label=="Anseranas_semipalmata")) #placed between Cygnus and Branta, node 7
# 
# #rename tips to genus only
# tree2$tip.label<-sapply(strsplit(tree2$tip.label,"_"),function(x) x[1])
# write.tree(tree2)
# 
# #family only tree
# tree <- read.newick("jetztreecleaned.newick")
# tips <- tree$tip.label
# 
# jetzspgenfam <- read_delim("//franklin.anu.edu.au/Home/u6562250/Teaching/BIOL3208/Maisie/3208 Final codes and Spreadsheets/Analysis/jetzspgenfam.csv", 
#                            "\t", escape_double = FALSE, trim_ws = TRUE)
# colnames(jetzspgenfam) <- c("PHYLO_NAME","JETZ_GEN","FAMILY")
# jetzspgenfam$PHYLO_NAME <- gsub(" ", "_", jetzspgenfam$PHYLO_NAME)
# 
# gentree <- read.newick("jetzgenSPtreecleaned.newick") #genus tree but tip names are species
# 
# jetzfam <- subset(jetzspgenfam, jetzspgenfam$PHYLO_NAME%in%gentree$tip.label) #only species that were chosen in the genus tree
# 
# unique(jetzfam$FAMILY)
# 
# jetzfam1 <- jetzfam[!duplicated(jetzfam$FAMILY),] #no duplicates, one from each family
# 
# #tree pruned to only selected species
# tips <- data.frame(name=tree$tip.label) #tips in original tree
# tips1 <- subset(tips, !tips$name%in%jetzfam1$PHYLO_NAME) #tips that we don't want 
# #RECHECK We are missing 326
# drop <- as.character(droplevels(tips1$name))
# 
# famtree <- drop.tip(tree, drop) #tree with only what we want
# unique(famtree$tip.label) #one species from each family
# #write.tree(famtree)
# 
# #renaming tips to family names
# match <-  data.frame(tips=famtree$tip.label)
# match <- merge(match, jetzfam1, by.x="tips", by.y="PHYLO_NAME")
# match <- match[,c(1,3)]
# 
# famtree1 <- sub.taxa.label(famtree, match)
# famtree1$tip.label

#BETADIV CALCULATIONS
###########PHYLOGENOMIC BETA DIVERSITY#######
#species
pair.j <- beta.pair(iaabird.sp, index.family="jaccard")
pair.j

#####Genus Beta Diversity####
pair.g.j <- beta.pair(iaabird.gen, index.family="jaccard")

######Family Beta Diversity######
pair.f.j <- beta.pair(iaabird.f, index.family="jaccard")

betadiv.sum <- rbind(beta.multi(iaabird.sp, index.family = "jaccard"), 
                     beta.multi(iaabird.gen, index.family = "jaccard"), 
                     beta.multi(iaabird.f, index.family = "jaccard"))
betadiv.sum <- as.data.frame((betadiv.sum))
betadiv.sum$level <- c("Species", "Genus", "Family")

#dendrogram plotting
par(mar=c(10,5,2,2))
region = data.frame(number = 1:11)
region$region <- c("Sahul","Sahul","Island", "Island", "Island", "Island", "Island", "Sunda", "Sunda", "Sunda", "Sunda")
regiontype <- factor(region$region)
n_region_types <- length(unique(regiontype))
cols_6 <- c("#cccccc", "#fde725", "#55cc66")
col_region_type <- cols_6[regiontype]

# Now we can use them
hc.s <- hclust(pair.j$beta.jac)
# hc.s <- with(region, reorder(hc.s, number))
dend.s <- as.dendrogram(hc.s)
dend.s <- reorder(dend.s, region$number)
plot(dend.s, ylab = expression(beta[jac]))
abline(a = mean(pair.j$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")
# plot(hc.s, ylab = expression(beta[jac]), hang=-1)

k234 <- cutree(hc.s, k = 2:6)

title(main = "a) Species")
colored_bars(cbind(k234[,5:1], col_region_type), dend.s, rowLabels = c(paste0("k = ", 6:2), "Region"))

hc.g <- hclust(pair.g.j$beta.jac)
#hc.g <- with(region, reorder(hc.g, number))
dend.g <- as.dendrogram(hc.g)
dend.g <- reorder(dend.g, region$number)
plot(dend.g, ylab = expression(beta[jac]))

k234 <- cutree(hc.g, k = 2:6)
title(main = "b) Genus")
colored_bars(cbind(k234[,5:1], col_region_type), dend.g, rowLabels = c(paste0("k = ", 6:2), "Region"))

hc.f <- hclust(pair.f.j$beta.jac)
#hc.f <- with(region, reorder(hc.f, number))
dend.f <- as.dendrogram(hc.f)
plot(dend.f, ylab = expression(beta[jac]))
k234 <- cutree(hc.f, k = 2:6)
title(main = "c) Family")
colored_bars(cbind(k234[,5:1], col_region_type), dend.f, rowLabels = c(paste0("k = ", 6:2), "Region"))

par(mfrow=c(3,1))
par(mar=c(5,2,2,2))
plot(dend.s, xlab = expression(beta[jac]), horiz = T)
abline(v = mean(pair.j$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")
title(main = "a) Species")

plot(dend.g, xlab = expression(beta[jac]), horiz = T)
abline(v = mean(pair.g.j$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")
title(main = "b) Genus")

plot(dend.f, xlab = expression(beta[jac]), horiz = T)
abline(v = mean(pair.f.j$beta.jac), b = 0, lwd = 3, lty = 2, col = "Red")
title(main = "c) Family")

#as matrix for each beta diversity index
spDist <-  as.matrix(pair.j$beta.jac)
gDist <-  as.matrix(pair.g.j$beta.jac)
fDist <-as.matrix(pair.f.j$beta.jac)

#VISUALIZING BETADIV MATRIX
spDist.re <- melt(round(spDist, digits =2))
sp.heat <- ggplot(spDist.re, aes(Var1, Var2, fill=value)) + geom_tile() + 
        geom_text(aes(label = value)) + 
        scale_fill_viridis(name= expression(paste(beta, "-diversity")), breaks=c(0,0.5,1),labels=c(0,0.5,1), limits=c(0,1)) +
        scale_x_discrete(limits=c("AUS","WMN","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA"))+
        scale_y_discrete(limits=c("AUS","WMN","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA"))+
        xlab(" ") + ylab(" ") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
sp.heat

gDist.re <- melt(round(gDist, digits =2))
g.heat <- ggplot(gDist.re, aes(Var1, Var2, fill=value)) + geom_tile() + 
        geom_text(aes(label = value)) + 
        scale_fill_viridis(name= expression(paste(beta, "-diversity")), breaks=c(0,0.5,1),labels=c(0,0.5,1), limits=c(0,1)) +
        scale_x_discrete(limits=c("AUS","WMN","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA"))+
        scale_y_discrete(limits=c("AUS","WMN","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA"))+
        xlab(" ") + ylab(" ") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
g.heat

fDist.re <- melt(round(fDist, digits =2))
f.heat <- ggplot(fDist.re, aes(Var1, Var2, fill=value)) + geom_tile() + 
        geom_text(aes(label = value)) + 
        scale_fill_viridis(name= expression(paste(beta, "-diversity")), breaks=c(0,0.5,1),labels=c(0,0.5,1), limits=c(0,1)) +
        scale_x_discrete(limits=c("AUS","WMN","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA"))+
        scale_y_discrete(limits=c("AUS","WMN","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA"))+
        xlab(" ") + ylab(" ") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
f.heat

#VISUALIZING CONNECTIVITY ANALYSIS
nodes <- c("AUS","NG", "EMN", "MAU", "SUL", "PHI","LSU","BOR","JAV","SUM","SEA")
nodes <- as.data.frame(nodes)
nodes$id <-  c(1:11)
nodes$sp.r <- c(753,740, 499, 406, 408, 577, 444, 581, 477, 619, 879) #species
nodes$gen.r <- c(338, 318, 188, 197, 211, 254, 209, 298, 270, 310, 375) #genus
nodes$fam.r <- c(95, 92, 70, 74, 75, 84, 77, 88, 92, 92, 95) #family
nodes$shelf <- c(rep(1, 2), rep(2,5), rep(3,4))
nodes <-  nodes[,c(2,1,6, 3,4,5 )]
colnames(nodes) <- c("id", "country", "shelf","sp.r", "gen.r", "fam.r")
nodes$country <-  as.character(nodes$country)

shelf.col <- c(rep("#fde725",2), rep("#cccccc",5), rep("#55cc66",4))

#SPECIES CONNECTIVITY
edges.sp <- read_csv("edges.sp.csv")

colnames(edges.sp) <- c("Var1","Var2","weight")

net.sp  <- graph_from_data_frame(d=edges.sp, vertices = nodes, directed=T)

E(net.sp)$width <-  E(net.sp)$weight/8
E(net.sp)$arrow.mode <- 0
V(net.sp)$label <- c("AUS","NG","EMN", "MAU", "SUL", "PHI","LSU","BOR","JAV","SUM","SEA")
V(net.sp)$size <- V(net.sp)$sp.r/15
V(net.sp)$frame.color <- c(rep("#fde725",2), rep("#cccccc",5), rep("#55cc66",4))

weight.sum <-  sum(edges.sp$weight)
cut.low <- 0.01*weight.sum
cut.mid <- 0.02*weight.sum

E(net.sp)$color <- ifelse( E(net.sp)$weight <= cut.low, "#dcdcdc","#4d4d4d")

l <-  layout_in_circle(net.sp.n)
plot(net.sp, layout=l, 
     #edge.color="gray85",
     vertex.label.font= 2,
     vertex.label.color="black",
     vertex.label.cex=.7,
     vertex.color=shelf.col)

#GENUS CONNECTIVITY
edges.g <- read_csv("edges.g.csv")

colnames(edges.g) <- c("Var1","Var2","weight")

net.g  <- graph_from_data_frame(d=edges.g, vertices = nodes, directed=T)

E(net.g)$width <-  E(net.g)$weight/5
E(net.g)$arrow.mode <- 0
V(net.g)$label <- c("AUS","NG","EMN", "MAU", "SUL", "PHI","LSU","BOR","JAV","SUM","SEA")
V(net.g)$size <- V(net.g)$gen.r/7
V(net.g)$frame.color <- c(rep("#fde725",2), rep("#cccccc",5), rep("#55cc66",4))

weight.sum <-  sum(edges.g$weight)
cut.low <- 0.01*weight.sum

E(net.g)$color <- ifelse( E(net.g)$weight <= cut.low, "#dcdcdc","#4d4d4d")

l <-  layout_in_circle(net.g.n)
plot(net.g, layout=l, 
     #edge.color="gray85",
     vertex.label.font= 2,
     vertex.label.color="black",
     vertex.label.cex=.7,
     vertex.color=shelf.col)

#FAMILYCONNECTIVITY
edges.f <- read_csv("edges.f.csv")

colnames(edges.f) <- c("Var1","Var2","weight")

net.f  <- graph_from_data_frame(d=edges.f, vertices = nodes, directed=T)

E(net.f)$width <-  E(net.f)$weight
E(net.f)$arrow.mode <- 0
V(net.f)$label <- c("AUS","NG","EMN", "MAU", "SUL", "PHI","LSU","BOR","JAV","SUM","SEA")
V(net.f)$size <- V(net.f)$fam.r/2
V(net.f)$frame.color <- c(rep("#fde725",2), rep("#cccccc",5), rep("#55cc66",4))

weight.sum <-  sum(edges.f$weight)
cut.low <- 0.01*weight.sum

E(net.f)$color <- ifelse( E(net.f)$weight <= cut.low, "#dcdcdc","#4d4d4d")

l <-  layout_in_circle(net.f.n)
plot(net.f, layout=l, 
     #edge.color="gray85",
     vertex.label.font= 2,
     vertex.label.color="black",
     vertex.label.cex=.7,
     vertex.color=shelf.col)

####### GEOG ENV DISTANCe
env <- read_csv("biomeans.csv")
#determine best set of non correlated variable
thr <- 0.8 
n <- 5
r <- cor(env)
score <- 1:10
combos <- combn(ncol(r), n)
summed_score <- apply(combos, 2, function(x) {
        z <- abs(r[x, x])
        if(any(z[lower.tri(z)] > thr)) NA else sum(score[x])
})

min(summed_score, na.rm=T)
## [1] 13

which.min(summed_score)
## [1] 9

r[combos[,229], combos[,229]]
corrplot(r[combos[,229], combos[,229]], method = "number")

env <- env[,c("bio10","bio1","bio13","bio14","bio15")]

EnvDist <- as.matrix(dist(env, method="euclidean"))
EnvDist <- EnvDist/max(EnvDist)

#shelf matrix
ShelfCat <- read_csv("shelfmatrix.csv")
ShelfCat <- ShelfCat[,-1]
ShelfCat <- as.matrix(ShelfCat)
rownames(ShelfCat) <- c("AUS","NWG","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA")

#minimal distance matrix
m.GeogDist <- read.csv("region.shp.distance.csv")
m.GeogDist <- m.GeogDist[,-1]
m.GeogDist <- as.matrix(m.GeogDist)
m.GeogDist <- m.GeogDist/max(m.GeogDist)

Area.mx <- read.csv("area.mx.csv")
Area.mx <- Area.mx[,-1]
Area.mx <- Area.mx/max(Area.mx)
Area.mx <- as.matrix(Area.mx)
rownames(Area.mx) <- c("AUS","NWG","EMN","MAU","SUL","PHI","LSU","BOR","JAV","SUM","SEA")


#####
# MMRR performs Multiple Matrix Regression with Randomization analysis
# Y is a dependent distance matrix
# X is a list of independent distance matrices (with optional names)

MMRR<-function(Y,X,nperm=999){
        #compute regression coefficients and test statistics
        nrowsY<-nrow(Y)
        y<-unfold(Y)
        if(is.null(names(X)))names(X)<-paste("X",1:length(X),sep="")
        Xmats<-sapply(X,unfold)
        fit<-lm(y~Xmats)
        coeffs<-fit$coefficients
        summ<-summary(fit)
        r.squared<-summ$r.squared
        tstat<-summ$coefficients[,"t value"]
        Fstat<-summ$fstatistic[1]
        tprob<-rep(1,length(tstat))
        Fprob<-1
        
        #perform permutations
        for(i in 1:nperm){
                rand<-sample(1:nrowsY)
                Yperm<-Y[rand,rand]
                yperm<-unfold(Yperm)
                fit<-lm(yperm~Xmats)
                summ<-summary(fit)
                Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
                tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
        }
        
        #return values
        tp<-tprob/(nperm+1)
        Fp<-Fprob/(nperm+1)
        names(r.squared)<-"r.squared"
        names(coeffs)<-c("Intercept",names(X))
        names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
        names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
        names(Fstat)<-"F-statistic"
        names(Fp)<-"F p-value"
        return(list(r.squared=r.squared,
                    coefficients=coeffs,
                    tstatistic=tstat,
                    tpvalue=tp,
                    Fstatistic=Fstat,
                    Fpvalue=Fp))
}

# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR

unfold<-function(X){
        x<-vector()
        for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
        return(x)
}
####
l.matrix <- list(m.GeogDist, EnvDist, Area.mx, ShelfCat)

mmrr.df <- MMRR(spDist, l.matrix)
mmrr.df
mmrr.df2 <- MMRR(gDist, l.matrix)
mmrr.df2
mmrr.df3 <- MMRR(fDist, l.matrix)
mmrr.df3


all.data <- cbind(melt(spDist), melt(gDist), melt(fDist), 
                  melt(m.GeogDist), melt(EnvDist), melt(Area.mx), melt(ShelfCat))
all.data <- all.data[,c(1,2,3,6,9, 12, 15, 18, 21)]
colnames(all.data) <- c("Var1", "Var2", "spDist", "gDist", "fDist", "m.GeogDist", "EnvDist", "Area.mx", "Shelf")
write.csv(all.data, "all.dist.matrix.melted.csv")

all.data <- subset(all.data, all.data$spDist!=0)

all.data$Shelf <- as.factor(all.data$Shelf)

all.data$out <- "no"
all.data$out[all.data$Var1 == "AUS" & all.data$Var2 == "NWG"] <- "yes"
all.data$out[all.data$Var1 == "NWG" & all.data$Var2 == "AUS"] <- "yes"

sp.mGeog <- 
        ggplot(all.data, aes(x = m.GeogDist, y = spDist, color = Shelf)) +
        geom_point(size = 3) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        geom_abline(intercept = mmrr.df$coefficients[1], slope = mmrr.df$coefficients[2]) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab(" ")+
        ylab("Dissimilarity") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle("a) Species")
        
sp.Env <- 
        ggplot(all.data, aes(EnvDist, spDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df$coefficients[1], slope = mmrr.df$coefficients[3]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab(" ")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle(" ")

sp.Area <- 
        ggplot(all.data, aes(Area.mx, spDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df$coefficients[1], slope = mmrr.df$coefficients[4]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab(" ")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none")+
        ggtitle(" ")

sp.Shelf <- 
        ggplot(all.data, aes(Shelf, spDist, fill = Shelf)) +
        geom_boxplot() +
        geom_point(stat="unique", position = "jitter", size = 3) +
        scale_fill_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_discrete(labels = c("No", "Yes")) +
        xlab(" ")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none")+
        ggtitle(" ")

g.mGeog <- 
        ggplot(all.data, aes(m.GeogDist, gDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df2$coefficients[1], slope = mmrr.df2$coefficients[2]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab(" ")+
        ylab("Dissimilarity") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle("b) Genus")

g.Env <- 
        ggplot(all.data, aes(EnvDist, gDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df2$coefficients[1], slope = mmrr.df2$coefficients[3]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab(" ")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle(" ")

g.Area <- 
        ggplot(all.data, aes(Area.mx, gDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df2$coefficients[1], slope = mmrr.df2$coefficients[4]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab(" ")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none")+
        ggtitle(" ")

g.Shelf <- 
        ggplot(all.data, aes(Shelf, gDist, fill = Shelf)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(stat="unique", position = "jitter", size = 3) +
        scale_fill_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_discrete(labels = c("No", "Yes")) +
        xlab(" ")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none")+
        ggtitle(" ")

f.mGeog <- 
        ggplot(all.data, aes(m.GeogDist, fDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df3$coefficients[1], slope = mmrr.df3$coefficients[2]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab("Minimum Geographic Distance")+
        ylab("Dissimilarity") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle("c) Family")

f.Env <- 
        ggplot(all.data, aes(EnvDist, fDist, color = Shelf)) +
        geom_point(size = 3) +
        geom_abline(intercept = mmrr.df3$coefficients[1], slope = mmrr.df3$coefficients[3]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab("Environmental Distance")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle(" ")

f.Area <- 
        ggplot(all.data, aes(Area.mx, gDist, color = Shelf)) +
        geom_point(size = 3)+
        geom_abline(intercept = mmrr.df3$coefficients[1], slope = mmrr.df3$coefficients[4]) +
        scale_color_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_continuous(limits=c(0,1), expand = c(0,0)) +
        xlab("Total Land Area")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none")+
        ggtitle(" ")

f.Shelf <- 
        ggplot(all.data, aes(Shelf, fDist, fill = Shelf)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(stat="unique", position = "jitter", inherit.aes = TRUE, size = 3) +
        scale_fill_manual(values = c("#545454","#55cc66")) +
        scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
        scale_x_discrete(labels = c("No", "Yes")) +
        xlab("Connected during the last glacial maximum")+
        ylab(" ") +
        theme_minimal(base_size = 14)+
        theme(axis.line = element_line(colour = "black"), legend.position = "none")+
        ggtitle(" ")

grid.arrange(sp.mGeog, sp.Env, sp.Area, sp.Shelf,
             g.mGeog, g.Env, g.Area, g.Shelf, 
             f.mGeog, f.Env, f.Area, f.Shelf, nrow = 3)

sparea <- read.csv("./sparea.csv")
colnames(sparea) <- c("area", "size", "endemic_sp", "sp", "gen", "fam", "group")

ggplot(sparea, aes(x=log(size), y=(sp), label = area)) +
        geom_smooth(method = lm) +
        geom_point() +
        geom_text(hjust=-0.2, vjust=0) +
        theme_minimal() +
        xlab(expression(paste("log"[10]," (Area in km"^2,")"))) +
        ylab("Species Richness")


sparea$area <- factor(sparea$area, levels = sparea$area[order(sparea$size)])

sparea.size <- 
        ggplot(sparea, aes(x=area, y=log(size), fill=group)) +
        geom_bar(stat="identity") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("#cccccc", "#fde725", "#55cc66")) +
        xlab("Area") +
        ylab(expression(paste("log"[10]," (Area in km"^2,")"))) +
        theme_minimal() +
        theme(legend.position = "none")

sparea.sp <- 
        ggplot(sparea, aes(x=area, y=(sp), fill=group)) +
        geom_bar(stat="identity") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("#cccccc", "#fde725", "#55cc66")) +
        xlab("Area") +
        ylab("Species Richness") +
        theme_minimal() +
        theme(legend.position = "none") 

grid.arrange(sparea.size, sparea.sp, nrow = 1)

sparea.sum <- aggregate(sparea[, 2:5], list(sparea$group), mean)
sparea[sparea$area=="EMN",2:6] + sparea[sparea$area=="NG",2:6]

sparea.new <- rbind(sparea, c("NG", 761830, 637, 1241, 506, 102, "sahul"))
sparea.new <- sparea.new[-11,]
sparea.new <- sparea.new[-5,]
sparea.new$logsize <- log(as.numeric(sparea.new$size))

sparea.size.n <- 
        ggplot(sparea.new, aes(x=area, y=logsize, fill=group)) +
        geom_bar(stat="identity") +
        scale_y_continuous(expand=c(0,0)) +
        scale_fill_manual(values=c("#cccccc", "#fde725", "#55cc66")) +
        xlab("Area") +
        ylab(expression(paste("log"[10]," (Area in km"^2,")"))) +
        theme_minimal() +
        theme(legend.position = "none")
sparea.sp.n <- 
        ggplot(sparea.new, aes(x=area, y=as.numeric(sp), fill=group)) +
                geom_bar(stat="identity") +
                scale_y_continuous(expand=c(0,0)) +
                scale_fill_manual(values=c("#cccccc", "#fde725", "#55cc66")) +
                xlab("Area") +
                ylab("Species Richness") +
                theme_minimal() +
                theme(legend.position = "none") 
grid.arrange(sparea.size.n, sparea.sp.n, nrow = 1)

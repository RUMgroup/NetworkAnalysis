#Network Analysis: Bridging R and Cytoscape

#Aims
#Introduce igraph for basic network constuction and analysis
#Demonstrate the use of RCy3 and the REST API to visualise igraph networks in Cytoscape


####Intro####
#Why use networks?

#Allow modelling of interactions between biological entities
#Framework to intergrate information 
#-can add expression data onto a protein interaction network to find active pathways in disease
#-if you know some causative disease genes can you infer others by network properties?


#How do we do it in R?
#- Plethora of packages on CRAN/Biocondctor for network analysis
#igraph is popular and performs most basic graph functions, plus have many algorithms pre-implemented 
#also available for python


###igraph#####

library("igraph")
library("RCy3")
library("gplots")

#First make an edge list
g <- as.data.frame(matrix(c(1, 2, 1, 3, 1, 4, 3, 4, 2, 4),byrow = TRUE, ncol = 2))

#make the graph
g<-graph.data.frame(g,directed = F)

#graph.adjacency() can be useful for input as well

#inspect the graph object
g
summary(g)

#visualise - don't do for big networks!
plot(g)

#add a node 
g<-add.vertices(graph = g,nv = 1,name=5) #note we are returning the modified graph

#add an edge from 1 to 5
g<-add.edges(graph = g,edges = c(1,5,5,1))

plot(g)

#let's pipe 'cause we can
g<-as.data.frame(matrix(c(1, 2, 1, 3, 1, 4, 3, 4, 2, 4),byrow = TRUE, ncol = 2)) %>%
  graph.data.frame(directed = T) %>%
  add.vertices(1,name=5) %>%
  add.edges(c(1,5,5,1,4,3)) %>%
  set_edge_attr("color", value = "red")

plot(g)

#iterate over nodes and edges
V(g)
E(g)

#Can create new vertice and edge attributes 
V(g)$size<-50
V(g)$color<-"blue"

#the plot function will interoret certain attributes by default
plot(g)

#can assign attribute values to indivdual nodes
V(g)$color<-c("blue","red","orange","purple","green")
plot(g)

#can do more advanced visualisation
plot(g,mark.groups = list(c(1,5),c(2,3,4)))
     
#basic functions provide building blocks for advanced custom algorithms
#i.e
neighbors(g,v = 5,mode = "all")

#best to avoid numbers as node names, things get confusing with index vs names for functions

###### Sending graphs to cytoscape #####

#Method 1: - pretty manual
#need some function utilities - borrowed from idekar lab (https://github.com/idekerlab/cy-rest-R)
source("cytoscape_util.R")
source("utils.R")

#convert graph to json
g.json<-toCytoscape(g)

network.url = paste(base.url, "networks", sep="/")
res <- POST(url=network.url, body=g.json, encode="json")

# Extract network SUID from the return value
network.suid = unname(fromJSON(rawToChar(res$content)))

style.simple <- buildStyleSimple("simple", g.json, color = "color")
style.url = paste(base.url, "styles", sep="/")
POST(url=style.url, body=style.simple, encode = "json")

apply.style.url = paste(base.url, "apply/styles/simple", toString(network.suid), sep="/")
GET(apply.style.url)

layout.params = list(
  name="unweighted",
  value=TRUE
)
layout.params.url = paste(base.url, "apply/layouts/kamada-kawai/parameters", sep="/")
PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")

# Apply layout
params <- paste(toString(network.suid), "?column=id", sep="")
apply.layout.url = paste(base.url, "apply/layouts/kamada-kawai", params, sep="/")
GET(apply.layout.url)
                  

##Method 2 - less manual

#RCy3 only works with graphNEL objects so need to convert

g.NEL<-igraph.to.graphNEL(g)

g.NEL<-initNodeAttribute(g.NEL,"color","char","black")
g.NEL<-initNodeAttribute(g.NEL,"size","numeric",1)
g.NEL<-initEdgeAttribute(g.NEL,"color","char","black")
g.NEL<-initEdgeAttribute(g.NEL,"weight","char","black")



nodeDataDefaults(g.NEL,"label")<-"default"
nodeData(self = g.NEL,V(g)$name,"label")<-V(g)$name
g.NEL<-initNodeAttribute(g.NEL,"label","char",default.value = "name")

#send graph to cytoscape
g.Cyto<-CytoscapeWindow(title = "example",graph = g.NEL,overwriteWindow = T)
displayGraph(g.Cyto)

#set layout
layoutNetwork(g.Cyto,"force-directed")

#set color
setNodeColorDirect(g.Cyto,V(g)$name,V(g)$color)

#set node shape
setNodeShapeDirect(g.Cyto,V(g)$name,rep("triangle",5))

#The connection is two ways - can get info back from cytoscape
getSelectedNodes(g.Cyto)

#could run an app in cytoscape then do extra stats in R



#### Structural Network Analysis ####

#Let's look at the yeast protein interaction network

# Load yeast network SIF file as Data Frame
yeast.table <- read.table("yeastHighQuality.sif")

# Convert it to simple edge list
yeast.table.edgelist <- yeast.table[c(1,3)]

# Convert data frame to undirected igraph object
g.original <- graph.data.frame(yeast.table.edgelist, directed=F)

# Extract componentes (individual connected subgraphs)
subgraphs <- decompose.graph(g.original)

# Pick largest subgraph
largest.subgraph <- subgraphs[[which.max(sapply(subgraphs, vcount))]]

# Remove duplicate edges
g <- simplify(largest.subgraph, remove.multiple=T, remove.loops=T)


# Global Network Statistics
g$density <- graph.density(g) # Density
g$transitivity <- transitivity(g) # Transitivity

# Node statistics
V(g)$closeness <- closeness(g) # Closeness Centrarity
V(g)$degree <- igraph::degree(g) # Degree
V(g)$pagerank <- page.rank(g, directed = FALSE)$vector # PageRank
V(g)$betweenness <- betweenness(g) # Betweenness Centrarity

# Edge statistics
E(g)$betweenness.edge <- edge.betweenness(g) # Edge Betweenness

#Community Detection: Try multiple algorithms

#optimise objective function - modularity
communities.greedy <- fastgreedy.community(g)

#by eginvector of the modularity matrix
communities.leading <- leading.eigenvector.community(g)

#neighbour majority voting
communities.label.propagation <- label.propagation.community(g)

V(g)$community.greedy <- communities.greedy$membership
V(g)$community.leading <- communities.leading$membership
V(g)$community.label.propagation<- communities.label.propagation$membership

#add colour for the community labels

V(g)$colors.community.greedy <- communityToColors(
  communities.greedy$membership,
  length(communities.greedy))

V(g)$colors.community.leading <- communityToColors(
  communities.leading$membership,
  length(communities.leading))

V(g)$colors.community.label.propagation <- communityToColors(
  communities.label.propagation$membership,
  length(communities.label.propagation))

#work out which edges are intra community
E(g)$community.greedy <- getCommunityEdge(g, V(g)$community.greedy)
E(g)$community.leading <- getCommunityEdge(g, V(g)$community.leading)
E(g)$community.label.propagation <- getCommunityEdge(g, V(g)$community.label.propagation)

#color the edges by community
E(g)$colors.community.greedy <- communityToColors(array(E(g)$community.greedy), length(communities.greedy))
E(g)$colors.community.leading <- communityToColors(array(E(g)$community.leading), length(communities.leading))
E(g)$colors.community.label.propagation <- communityToColors(array(E(g)$community.label.propagation), length(communities.label.propagation))


#send it to cytoscape
g.NEL<-igraph.to.graphNEL(g)

#too many attributes to do one by one!
summary(g)

#wrote util functions to deal with it
g.NEL<-processEdgeAttributes(g,g.NEL)
g.NEL<-processNodeAttributes(g,g.NEL)


#send graph to cytoscape
g.Cyto<-CytoscapeWindow(title = "example",graph = g.NEL,overwriteWindow = T)
displayGraph(g.Cyto)

#Doesn't see to work and very slow node operations!

#set the background
setDefaultBackgroundColor(g.Cyto,"black")

#set edge width
setEdgeLineWidthDirect(g.Cyto,edges,new.value = 2)

#set color of nodes
setNodeColorDirect(g.Cyto,V(g)$name,V(g)$colors.community.greedy)
setEdgeColorDirect(g.Cyto,edges,E(g)$colors.community.greedy)

setEdgeOpacityRule(g.Cyto,edge.attribute.name =community)

#set Border width
setNodeBorderWidthDirect(g.Cyto,V(g)$name,new.sizes = rep(0,vcount(g)))


###PLAN B #####

#Seems far to slow for a large network, let's try a more direct approach

# Convert igraph object into Cytoscape.js JSON
cyjs <- toCytoscape(g)

# POST it to Cytoscape
network.url = paste(base.url, "networks", sep="/")
res <- POST(url=network.url, body=cyjs, encode="json")

# Extract network SUID from the return value
network.suid = unname(fromJSON(rawToChar(res$content)))

# Generate Visual Styles
style.greedy <- buildStyle("greedy", g, colors = "colors.community.greedy", community="community.greedy")
style.leading <- buildStyle("leading", g, colors = "colors.community.leading", community="community.leading")
style.label.propagation <- buildStyle("label.propagation", g,
                                      colors = "colors.community.label.propagation", community="community.label.propagation")

style.url = paste(base.url, "styles", sep="/")
POST(url=style.url, body=style.greedy, encode = "json")
POST(url=style.url, body=style.leading, encode = "json")
POST(url=style.url, body=style.label.propagation, encode = "json")

# Apply a Style
apply.style.url = paste(base.url, "apply/styles/greedy", toString(network.suid), sep="/")
GET(apply.style.url)

# Tweak Layout parameters
layout.params = list(
  name="unweighted",
  value=TRUE
)

layout.params.url = paste(base.url, "apply/layouts/kamada-kawai/parameters", sep="/")
PUT(layout.params.url, body=toJSON(list(layout.params)), encode = "json")

# Apply layout
params <- paste(toString(network.suid), "?column=community.greedy", sep="")
apply.layout.url = paste(base.url, "apply/layouts/kamada-kawai", params, sep="/")
GET(apply.layout.url)

# Perform Edge Bundling
apply.bundling.url = paste(base.url, "apply/edgebundling", toString(network.suid), sep="/")
GET(apply.bundling.url)

# Toggle graphics details
lod.url = paste(base.url, "ui/lod", sep="/")
PUT(lod.url)



#could go on to look at Gene ontology enrichement of clusters, TF binding to promoters of proteins,  -> biological significance of modules
#Interactions between modules?

#Meaningful visualisation on a large network is hard. More common to see smaller sub-networks of interest to study in detail.










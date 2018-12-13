##########		Redes de Cáncer Integradas		#########
# @lunysska
# @dario
# @miguel
# @ana
#
# Dec, 2018

#Modificar la ruta asignada para los archivos
setwd( "/Users/virg/Downloads" )
source( 'util.R' )
source( 'medidas_centralidad.R' )

# Instalar paquetes y cargar sus librerías
usa_paquete( "igraph", "tcltk" )

# Leer y asignar al archivo csv de interacciones
# de genese en cancer de cerebro
#cancer <- read.csv("string_interactions_pulmon.csv",
#                   header = TRUE, as.is = TRUE)
nodos <- read.csv( "nodos_red_cancer.csv", 
                  header=TRUE, 
                  as.is=T )

vinculos <- read.csv( "enlaces_red_cancer.csv",
                      header=TRUE,
                      as.is=T )

# Generar red de interaccion
red_cancer <- graph_from_data_frame( d=vinculos, 
                                    vertices=nodos, 
                                    directed=FALSE )
red_cancer
head( red_cancer )

#obtenemos el grado de los nodos
grado_de_nodos <- degree( red_cancer )

#Cambiamos el tamaño de los nodos, tomando en cuenta la proporción
#del grado máximo (74), y considerando el tamaño máximo que deseamos
#en este caso 20
V( red_cancer )$size <- (grado_de_nodos*20)/74


# Pimpeando el grafo
plot(red_cancer, 
     vertex.color="lightsteelblue2",
     vertex.frame.color="black",
     vertex.label.color = "black",
     edge.width=1, 
     edge.color= "gray85")

red_cancer
print_all(red_cancer)
class(red_cancer)
mode(red_cancer)
attributes(red_cancer)

# Abre una ventana que permite manipular el
# grafo
tkplot(red_cancer)

# Para hacer una matriz de la red de cancer
as_adj(red_cancer) # 157 x 157 sparse Matrix of class "dgCMatrix"
red_cancer_matriz <- as_adj(red_cancer)
red_cancer_matriz 

# dim -> para ver la dimension de la tabla
dim(red_cancer_matriz) # 	157	157

E(red_cancer) # enlista los vinculos de la red	773/773 edges
V(red_cancer) # enlista los nodos de la red		157/157 vertices

# Obtener atributos de los vinculos: names 
# Por omision extrae el atributo “name” del “from” “to” (primeras dos columnas) del archivo
V(red_cancer)$name	

# red simplificada completa "red_cancer_simpl" <- simplify
# No bucles, no direcciones y no multimvinculos
red_cancer_simple <- simplify(red_cancer)
red_cancer_simple
grado_de_nodos <- degree( red_cancer_simple )
V( red_cancer_simple )$size <- (grado_de_nodos*20)/74
plot( red_cancer_simple, 
     vertex.color="lightsteelblue2",
     vertex.frame.color="black",
     vertex.label.color = "black",
     edge.width=1, 
     edge.color= "gray85" )

# Para remover multivinculos y loops como si fuera directo "simplify"
red_cancer_simple <- simplify(red_cancer, remove.multiple = TRUE,
                              remove.loops = TRUE)
red_cancer_simple # IGRAPH c3c034b UN-- 157 769
head(red_cancer_simple)
grado_de_nodos <- degree( red_cancer_simple )
V( red_cancer_simple )$size <- (grado_de_nodos*20)/74
plot( red_cancer_simple, 
      vertex.color="lightsteelblue2",
      vertex.frame.color="black",
      vertex.label.color = "black",
      edge.width=1, 
      edge.color= "gray85" )

# Generar un grafo a partir de la matriz de adyacencia y la asigno como
# "red_cancer_grafo"
red_cancer_grafo <- graph_from_adjacency_matrix(red_cancer_matriz,
                                                mode = "directed",
                                                weighted = TRUE,
                                                diag = TRUE)
class(red_cancer_grafo) # "igraph"
red_cancer_grafo
# IGRAPH cc5a9ba DNW- 157 1538 -- 
# + attr: name (v/c), size_between (v/n), size_degree (v/n), size_eigen
# | (v/n), comp (v/n), weight (e/n)
# + edges from cc5a9ba (vertex names):
grado_de_nodos <- degree( red_cancer_grafo )
V( red_cancer_grafo )$size <- ( grado_de_nodos*20 )/74
plot( red_cancer_grafo, 
      vertex.color="lightsteelblue2",
      vertex.frame.color="black",
      vertex.label.color = "black",
      edge.width=1, 
      edge.color= "gray85" )

head(red_cancer_grafo) # 6 x 157 sparse Matrix of class "dgCMatrix"

# Obtener atributos de los vinculos
V(red_cancer_grafo) # + 157/157 vertices
V(red_cancer_grafo)$name
V(red_cancer_grafo)$size_between
V(red_cancer_grafo)$size_degree
V(red_cancer_grafo)$size_eigen
V(red_cancer_grafo)$comp
V(red_cancer_grafo)$weight

# Obtener atributos de los enlaces
E(red_cancer_grafo) # + 1538/1538 edges

# Genero un grafo no dirigido a partir de la matriz de adyacencia y la asigno
# como "red_cancer_grafo_undirected" no direccionada
red_cancer_grafo_undirected <- graph_from_adjacency_matrix(red_cancer_matriz,
                                                           mode = "undirected")
red_cancer_grafo_undirected
grado_de_nodos <- degree( red_cancer_grafo_undirected )
V( red_cancer_grafo_undirected )$size <- ( grado_de_nodos*20 )/74
plot( red_cancer_grafo_undirected, 
      vertex.color="lightsteelblue2",
      vertex.frame.color="black",
      vertex.label.color = "black",
      edge.width=1, 
      edge.color= "gray85" )


########### CENTRALIDADES
### Degree
grado_de_nodos <- degree( red_cancer_grafo )
sort( grado_de_nodos, decreasing = TRUE )
max( degree( red_cancer_grafo ) ) # 142
min( degree( red_cancer_grafo ) ) # 2
hist( degree( red_cancer_grafo ) )
hist( degree_distribution(red_cancer_grafo))
plot(degree_distribution(red_cancer_grafo))

### Betweenness
betweenness(red_cancer_grafo)
class( betweenness(red_cancer_grafo) )
max(betweenness(red_cancer_grafo)) # 7088.344
min(betweenness(red_cancer_grafo)) # 0

edge_betweenness(red_cancer_grafo)
max(edge_betweenness(red_cancer_grafo)) # 1698.005
min(edge_betweenness(red_cancer_grafo)) # 0.06267

### Eigenvector
eigen_centrality(red_cancer_grafo)
eigen_centrality(red_cancer_grafo)[1]
# Para obtener el nodo con mayor y menor centralidad de eigenvector
which.max(eigen_centrality(red_cancer_grafo)$vector) # TP53    5
which.min(eigen_centrality(red_cancer_grafo)$vector) # KRT15    74

### Closeness
closeness(red_cancer_grafo)
### Warning message:
### In closeness(red_cancer_grafo) :
### At centrality.c:2617 :closeness
### centrality is not well-defined for
### disconnected graphs
max(closeness(red_cancer_grafo)) # 0.001655629
min(closeness(red_cancer_grafo)) # 4.109139e-05

### Como se observa en el plot de (red_cancer_grafo) (descomentar para ver):
# plot(red_cancer_grafo)
### Un par de nodos se encuentran desvinculados del resto de los nodos:
### KRT13 -- KRT15
### Por lo tanto, se debe hacer un grafo del componente principal de la red, eliminando
### los nodos inconexos

### Para modificar el tamaño de los nodos de acuerdo a atributos
# Generar un atributo al nodo
# V(net)$size <- betweeness(net)  añadiendo "*" o "/" y un valor, p.e. 5

#V(red_cancer_grafo)$size_between <- betweenness(red_cancer_grafo)*4
# Al atributo "size" de los vertices le voy a asignar 
#la multiplicacion por 4 de red_cancer_grafo
#plot(red_cancer_grafo)

# Medida del grado (degree)
#(grafo_degree <- degree(red_cancer_grafo))

#V(red_cancer_grafo)$size_degree <- grafo_degree
# Al atributo "size" de los vertices le voy a asignar la multiplicacion 
#por 4 de gmed
#plot(red_cancer_grafo)

#red_degree <- plot(red_cancer_grafo, vertex.color="orangered2",
#                   vertex.frame.color="black", edge.color="gray",
#                   edge.width=4)

### Intermediacion (betweeness)
#(grafo_between <- betweenness(red_cancer_grafo, directed = FALSE))  

#V(red_cancer_grafo)$size_between <- grafo_between/3 # Dividir tamaÃ±o de nodos entre 3 porque estan alrededor de 15
#V(red_cancer_grafo)$size_between
#V(red_cancer_grafo)$name

#red_between <- plot(red_politica_grafo, vertex.color="orangered2",
#                    vertex.frame.color="black", edge.color="gray",
#                    edge.width=4)

### Cercania
#(grafo_close <- closeness(red_cancer_grafo))
### Warning message:
### In closeness(red_cancer_grafo) :
### At centrality.c:2617 :closeness centrality is not
### well-defined for disconnected graphs
#min(closeness(red_cancer_grafo)) # 4.109139e-05
#max(closeness(red_cancer_grafo)) # 0.001655629

#V(red_cancer_grafo)$size_close <- grafo_close*1000 # Por 1000 para poder ajustar el tamaÃ ±o, ya que su cercania tiene valores entre 0.04 y 0.02
#V(red_cancer_grafo)$size_close
#V(red_cancer_grafo)$name

#red_close <- plot(red_cancer_grafo, vertex.color="orangered2",
#                  vertex.frame.color="black", edge.color="gray",
#                  edge.width=4)


#############################
######  Cliques##############
#############################
#Una de las formas de poder encontrar comunidades en una red buscando 
#cuales de los nodos tienen mayor cantidad de vinculos para ello ocupamos 
#la funcion de cliques.
#Los cliques son subgraficos maximales. Cada vertice esta conectado a 
#todos los demas.
#Se genera una tabla con el numero de subgraficas que se generan con un 
#numero determinado de nodos
table(sapply(cliques(red_cancer_grafo),length))

#esta función genera el numero maximal de cliques es decir que uno no este contenido en
#el siguiente.
table(sapply(maximal.cliques(red_cancer_grafo), length))

###################################
##### Clustering para red_cancer_grafo
###################################
###################################
###################################
###################################
grafo_cancer2 <- as.undirected(red_cancer_grafo, mode = "collapse")
grafo_cancer2 # IGRAPH a5ee028 UNW- 157 769
plot(grafo_cancer2)

grafo_cancer3 <- as.undirected(red_cancer_grafo, mode = "mutual")
grafo_cancer3 # IGRAPH b9282a9 UNW- 157 769
plot(grafo_cancer2, vertex.size= 15)


##### Comunidades grafo_cancer2
### Algoritmos para etiquetar

### Modularity
# Valor numerico que calcula bajo la clasificacion que le dimos,
# el num de vinculos respecto a los vinculos con los otros 
# grupos -> modularidad alta
# Entre - 1/2 y 1

##### cluster_walktrap() para red_cancer_grafo
cw <- cluster_walktrap(grafo_cancer2)
modularity(cw) # [1] 0.3571618
# Ve cuantos vinculos hay entre las etiquetas que di, si te acercas a 1 la
# modularidad es alta y refleja la separacion de los grupos
# Si te acercas a 1 la forma en la que los etiquetaste es muy buena

membership(cw)
plot(cw, grafo_cancer2)

##### cluster_edge_betweenness() para red_cancer_grafo
ceb <- cluster_edge_betweenness(grafo_cancer2)
### Warning messages:
### 1: In cluster_edge_betweenness(gpol2) :
### At community.c:460 :Membership vector will be selected based
### on the lowest modularity score.
### 2: In cluster_edge_betweenness(gpol2) :
### At community.c:467 :Modularity calculation with weighted
### edge betweenness community detection might not make
### sense -- modularity treats edge weights as similarities
### while edge betwenness treats them as distances
modularity(ceb) # [1] 0.355115
membership(ceb)
plot(ceb, grafo_cancer2) # Les pone el mismo color y los agrupa

##### cluster_spinglass() para red_cancer_grafo
cs <- cluster_spinglass(grafo_cancer2)
### Warning messages:
### 1: In cluster_edge_betweenness(gpol2) :
### At community.c:460 :Membership vector will be selected based
### on the lowest modularity score.
### 2: In cluster_edge_betweenness(gpol2) :
### At community.c:467 :Modularity calculation with weighted
### edge betweenness community detection might not make
### sense -- modularity treats edge weights as similarities
### while edge betwenness treats them as distances
modularity(cs)
membership(cs)
plot(cs, grafo_cancer2)


##### cluster_fast_greedy() para red_cancer_grafo
cfg <- cluster_fast_greedy(grafo_cancer2)
modularity(cfg) # [1] 0.4180651
membership(cfg)
plot(cfg, grafo_cancer2)

##### cluster_label_prop() para red_cancer_grafo
clp <- cluster_label_prop(grafo_cancer2)
modularity(clp) # [1] 0.1444039
membership(clp)
plot(clp, grafo_cancer2)


##### Creando data frame para la red
# Extraer tabla de vertices y otra de vinculos en dos archivos distintos

nodos_grafo_cancer2 <- as_data_frame(grafo_cancer2, what = "vertices")
enlaces_grafo_cancer2 <- as_data_frame(grafo_cancer2, what = "edges")

write.csv(nodos_grafo_cancer2, file = "nodos_grafo_cancer2.csv", row.names = FALSE)
write.csv(enlaces_grafo_cancer2, file = "enlaces_grafo_cancer2.csv", row.names = FALSE)

nodes <- read.csv("nodos_red_cancer.csv", header = TRUE, as.is = TRUE)
links <- read.csv("enlaces_red_cancer.csv", header = TRUE, as.is = TRUE)

head(nodes) # name size_between size_degree size_eigen comp size_close
nodes$id <- nodes$name
head(nodes)

head(links) # from     to weight
links$id <- links$weight
head(links) # from     to weight id




##### Animacion de la red en Java

install.packages("visNetwork")
library(visNetwork) # Hace animaciones en Java

visNetwork(nodes, links, width = "100%", height = "400px") 
# Corre pero a mi no me aparece nada

visNetwork(nodes, links, width = "100%", height = "400px",
           background = "#eeefff", main = "Red", submain = "",
           footer = "Hyperlinks and mentions among media sources")
# Corre pero solo me aparece el titulo “Red” y el fondo azul claro


##############################
### Como aparece una advertencia en el anaisis de closeness en red_cancer_grafo, 
### dado que la centralidad no esta bien definida para grafos disconexos:


##############################

### Eliminar clusters aislados de red_cancer_grafo

##############################


cl <- clusters(red_cancer_grafo)
cl$membership # Vector donde el elemento 1 contiene el indice
# de cluster del nodo i-1
cl$csize # Vector donde el elemento i da el tamaño del cluster i # [1] 155   2
cl$no # Para obtener el numero de clusters # [1] 2

which(red_cancer_grafo$csize <= 3) # Para recuperar los indices de
# los clusters que tendremos que
# borrar # integer(0)

small_clusters <- which(cl$csize <= 2)
vertices_to_delete <- which(cl$membership %in% small_clusters) -1

grafo_conectado <- delete.vertices(red_cancer_grafo, vertices_to_delete)
plot(grafo_conectado) 

##### Esta estrategia no funcionó para eliminar los disconexos, por lo tanto:
##############################

### Obtener el componente principal de la red_cancer_grafo

##############################


# Obtener el componente principal de la red (main cluster)
V(red_cancer_grafo)$comp <- components(red_cancer_grafo)$membership

# Generar un subconjunto basado en esos componentes
main <- induced_subgraph(red_cancer_grafo, V(red_cancer_grafo)$comp==1)
plot(main)

main # IGRAPH f4d4561 DNW- 155 1536


# Para ver opciones de layout
?igraph::layout


#########################################

### Medidas de centralidad solo del main cluster de la red_cancer_grafo

### Degree solo del main cluster de la red_cancer_grafo
degree(main)
max(degree(main)) # [1] 142
min(degree(main)) # [1] 2
hist(degree(main))

hist(degree_distribution(main))
plot(degree_distribution(main))

### Betweenness solo del main cluster de la red_cancer_grafo
betweenness(main)
max(betweenness(main)) # [1] 7088.344
min(betweenness(main)) # [1] 0

edge_betweenness(main)
max(edge_betweenness(main)) # [1] 1698.005
min(edge_betweenness(main)) # [1] 0.06267806




### Eigenvector solo del main cluster de la red_cancer_grafo
eigen_centrality(main)
eigen_centrality(main)[1]

# Closeness solo del main cluster de la red_cancer_grafo
closeness(main)
max(closeness(main)) # [1] 0.003448276
min(closeness(main)) # [1] 0.0009049774

### Para modificar el tamaño de los nodos de acuerdo a atributos
# Generar un atributo al nodo
# V(net)$size <- betweeness(net) "*" o "/" algo

V(main)$size_between <- betweenness(main)*4
plot(main)

# Medida del grado (degree)
(main_grafo_degree <- degree(main))

V(main)$size_degree <- main_grafo_degree
# Al atributo "size" de los vertices le voy a asignar la multiplicacion por 4 de gmed
plot(main_grafo_degree)

main_red_degree <- plot(main, vertex.color="orangered2",
                        vertex.frame.color="black", edge.color="gray",
                        edge.width=0.05)


### Intermediacion (betweeness)
(main_grafo_between <- betweenness(main, directed = FALSE))

V(main)$size_between <- main_grafo_between/3 # Dividir tamaño de nodos entre 3 porque estan alrededor de 15
V(main)$size_between
V(main)$name

main_red_between <- plot(main, vertex.color="orangered2",
                         vertex.frame.color="black", edge.color="gray",
                         edge.width=4)

### Cercania
(main_grafo_close <- closeness(main))
min(closeness(main))
max(closeness(main))

V(main)$size_close <- main_grafo_close*1000 # Por 1000 para poder ajustar el tamaño, ya que su cercania tiene valores entre 0.04 y 0.02
V(main)$size_close
V(main)$name

main_red_close <- plot(main, vertex.color="orangered2",
                       vertex.frame.color="black", edge.color="gray",
                       edge.width=4)

### Centralidad de eigenvector
(main_grafo_eigen <- eigen_centrality(main))
(main_grafo_eigen <- eigen_centrality(main)[1])
class(main_grafo_eigen)

main_grafo_eigen$vector
class(main_grafo_eigen$vector)

main_grafo_eigen <- main_grafo_eigen$vector
main_grafo_eigen

main_grafo_eigen <- as.numeric(main_grafo_eigen)
main_grafo_eigen

V(main)$size_eigen <- main_grafo_eigen*15
V(main)$size_eigen

main_red_eigen <- plot(main, vertex.color="orangered2",
                       vertex.frame.color="black", edge.color="gray",
                       edge.width=4, vertex.label.color="black")


which.max(eigen_centrality(main)$vector)
V(main)$name[3]





###################################

##### cluster_walktrap()

###################################

main_grafo_cancer2 <- as.undirected(main, mode = "collapse")
main_grafo_cancer2
plot(main_grafo_cancer2)

main_grafo_cancer3 <- as.undirected(main, mode = "mutual")
main_grafo_cancer3
plot(main_grafo_cancer2, vertex.size= 15)

##### Comunidades grafo_cancer2
### Algoritmos para etiquetar

### Modularity
# Valor numerico que calcula bajo la clasificacion que le dimos,
# el num de vinculos respecto a los vinculos con los otros grupos -> modularidad alta
# Entre - 1/2 y 1


main_cw <- cluster_walktrap(main_grafo_cancer2)
modularity(main_cw)
# Ve cuantos vinculos hay entre las etiquetas que di, si te acercas a 1 la
# modularidad es alta y refleja la separacion de los grupos
# Si te acercas a 1 la forma en la que los etiquetaste es muy buena


membership(main_cw)

main_ceb <- cluster_edge_betweenness(main_grafo_cancer2)
### Warning messages:
### 1: In cluster_edge_betweenness(gpol2) :
### At community.c:460 :Membership vector will be selected based
### on the lowest modularity score.
### 2: In cluster_edge_betweenness(gpol2) :
### At community.c:467 :Modularity calculation with weighted
### edge betweenness community detection might not make
### sense -- modularity treats edge weights as similarities
### while edge betwenness treats them as distances
modularity(main_ceb)
membership(main_ceb)
plot(main_ceb, main_grafo_cancer2) # Les pone el mismo color y los agrupa

main_cs <- cluster_spinglass(main_grafo_cancer2)
modularity(main_cs)
membership(main_cs)
plot(main_cs, main_grafo_cancer2)

main_cfg <- cluster_fast_greedy(main_grafo_cancer2)
modularity(main_cfg)
membership(main_cfg)
plot(main_cfg, main_grafo_cancer2)

main_clp <- cluster_label_prop(main_grafo_cancer2)
modularity(main_clp)
membership(main_clp)
plot(main_clp, main_grafo_cancer2)


##### Creando data frame para la red
# Extraer tabla de vertices y otra de vinculos en dos archivos distintos

nodos_main_grafo_cancer2 <- as_data_frame(main_grafo_cancer2, what = "vertices")
enlaces_main_grafo_cancer2 <- as_data_frame(main_grafo_cancer2, what = "edges")

write.csv(nodos_main_grafo_cancer2, file = "nodos_main_grafo_cancer2.csv", row.names = FALSE)
write.csv(enlaces_main_grafo_cancer2, file = "enlaces_main_grafo_cancer2.csv", row.names = FALSE)

main_nodes <- read.csv("nodos_main_grafo_cancer2.csv", header = TRUE, as.is = TRUE)
main_links <- read.csv("enlaces_main_grafo_cancer2.csv", header = TRUE, as.is = TRUE)

head(main_nodes)
main_nodes$id <- main_nodes$name
head(main_nodes)

head(main_links)
main_links$id <- main_links$weight

############################################################################
#"k-core" 
############################################################################
# Estas son los paquetes que vamos a utilizar
install.packages("devtools")
install_github("r-lib/devtools")
library(devtools)
## A algunos les funciono los siguiente  #####
install_github("DougLuke/UserNetR")
library(UserNetR)
### Pero si no fue así prueba lo siguiente ###
options(unzip="internal")
devtools::install_github("DougLuke/UserNetR",force=TRUE)
library(UserNetR)

###################################################
cancer <- read.csv("Redes_Cancer.csv",
                   header = TRUE, as.is = TRUE)
cancer
head(cancer)

# Generar red de interaccion
red_cancer <- graph_from_data_frame(cancer, directed = FALSE)
red_cancer
plot(red_cancer)

######################################
# medimos la densidad de la red
graph.density(red_cancer) #0.063122

# Para identificar la estructura de "k-core" en la red, usamos la 
# funcion "graph.coreness( )". Esta regresa un vector que enlista los 
# "k-core" de la  red. 


coreness <- graph.coreness(red_cancer)
coreness

# table() genera una tabla de frecuencias
# El resultado nos dice el rango de "k-cores" de 1 a 13

table(coreness)
#coreness
# 1  2  3  4  5  6  8  9 10 11 12 13 
#23 15 20 24  7  5 12 19  4  5  4 19 

# Podemos ver que el 13-core esta formado de 
# 19 nodos, el 12-core de 4 nodos 

maxCoreness <- max(coreness)

maxCoreness #13 maxCorness

V(red_cancer)$name  
V(red_cancer)$color <- coreness +1
V(red_cancer)$color

# par(mfrow=c(1,2)) un renglon dos columnas, mar= C(abajo, izq, arrib, derecha)
# de abajo a favor de las manecilas del reloj

op <- par(mar=rep(0,4)) # sin marco
plot(red_cancer, vertex.label.cex=0.5)
par(op)

dev.off()

# usamos la paleta de colres rainbow()
# rainbow puede dar una paleta de colores con un solo 
# argumento que es el numero "n" de colores que se necesitan
# en nuestro caso generamos una paleta del mismo tamaño 
# que el maximo "coreness"

colors <- rainbow(maxCoreness)
colors # codigo RGB

op <- par(mar=rep(0,4)) # sin marco
plot(red_cancer, vertex.label=coreness, vertex.color=colors[coreness])
par (op)

dev.off()

# el k-core mas grande es el 4-core y esta compuesto de 24 de los
# 157 nodos totales.

# Vamos a ir gradualmente de los mas bajos a los más altos

V(red_cancer)$name <- coreness
V(red_cancer)$color <- colors[coreness]
red_cancer1_13 <- red_cancer
red_cancer2_13 <- induced.subgraph(red_cancer, vids=which(coreness>1))
red_cancer3_13 <- induced.subgraph(red_cancer, vids=which(coreness>2))
red_cancer4_13 <- induced.subgraph(red_cancer, vids=which(coreness>3))
red_cancer5_13 <- induced.subgraph(red_cancer, vids=which(coreness>4))
red_cancer6_13 <- induced.subgraph(red_cancer, vids=which(coreness>5))
red_cancer7_13 <- induced.subgraph(red_cancer, vids=which(coreness>6))
red_cancer8_13 <- induced.subgraph(red_cancer, vids=which(coreness>7))
red_cancer9_13 <- induced.subgraph(red_cancer, vids=which(coreness>8))
red_cancer10_13 <- induced.subgraph(red_cancer, vids=which(coreness>9))
red_cancer11_13 <- induced.subgraph(red_cancer, vids=which(coreness>10))
red_cancer12_13 <- induced.subgraph(red_cancer, vids=which(coreness>11))
red_cancer13_13 <- induced.subgraph(red_cancer, vids=which(coreness>12))

lay <- layout.fruchterman.reingold(red_cancer)
lay
op <- par(mfrow=c(4,2), mar=c(1,0,1,0))

plot(red_cancer1_13, layout=lay, main="Todos los k-cores")
plot(red_cancer2_13, layout=lay[which(coreness>1),], main="k-cores 2-13")
plot(red_cancer3_13, layout=lay[which(coreness>2),], main="k-cores 3-13")
plot(red_cancer4_13, layout=lay[which(coreness>3),], main="k-cores 4-13")
plot(red_cancer5_13, layout=lay[which(coreness>4),], main="k-cores 5-13")
plot(red_cancer6_13, layout=lay[which(coreness>5),], main="k-cores 6-13")
plot(red_cancer7_13, layout=lay[which(coreness>6),], main="k-cores 7-13")
plot(red_cancer8_13, layout=lay[which(coreness>7),], main="k-cores 8-13")
plot(red_cancer9_13, layout=lay[which(coreness>8),], main="k-cores 9-13")
plot(red_cancer10_13, layout=lay[which(coreness>9),], main="k-cores 10-13")
plot(red_cancer11_13, layout=lay[which(coreness>10),], main="k-cores 11-13")
plot(red_cancer12_13, layout=lay[which(coreness>11),], main="k-cores 12-13")
plot(red_cancer13_13, layout=lay[which(coreness>12),], main="k-cores 13-13")
par(op)

############################################################################
#    Finaliza ANA
############################################################################


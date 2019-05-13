
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1403924.svg)](https://doi.org/10.5281/zenodo.1403924)

To install, 

library(devtools)
install_github("abnormally-distributed/rsfcNet")

If you use this package please cite as follows: 

**APA Citation:**

*Vaughan, B. (2018, August 8). rsfcNet: An R Package for Resting State Functional Connectivity Network Analysis  Zenodo. http://doi.org/10.5281/zenodo.1403924*


<p><br />
<em><strong>Introduction</strong></em><br />
<br />
Graph theory is an important and versatile approach to studying the connectivity of different brain regions and can be applied to any number of neuroimaging methods, including structural MRI, functional MRI, diffusion tensor or kurtosis imaging, EEG, and MEG. A graph is simply a matrix that quantifies the connections between individual members of some group of interest. Each member is called a node (or vertex) and the connection between each node is an edge. Such graphs can be binary, where each edge is either a 1 or a 0, or weighted, where the strength of connection can vary. Furthermore, weighted networks can be signed, as in the case of a correlation matrix. Correlations (or partial correlations) quantify a primary type of connectivity between given pair of nodes, but when considering many variables higher-order relationships can illuminate aspects of the network in question that bi-variate relationships fail to capture. However, correlations are not the only means by which one can define edges in a network.</p>

<p>Most graph theoretic metrics are centrality measures, which aim to identify through various methods the most important members of a network. Which measure is most appropriate depends on the type of network being characterized, the modality of being important (for example, a CEO of a company may have fewer direct connections than others within a company but nevertheless occupy a privileged and important place in the company network), and the means by which edges are represented.&nbsp;</p>

<p><br />
<em><strong>Network Construction</strong></em><br />
<br />
<strong>Data Preparation</strong></p>

<p>The first step is importing the data and preparing it for analysis. Assuming one has both time series files and confound matrices from a preprocessing pipeline in FSL or some other software, the first step is to import both.&nbsp;<br />
First, create variables pointing to the files in the directory.&nbsp;</p>

<pre>
<code>ts_list = get_file_list("C:/Users/YourName/Documents/fMRIstudy/timeseries/")

confound_list = get_file_list("C:/Users/YourName/Documents/fMRIstudy/confounds/")</code></pre>

<p>Next import the files from the lists.</p>

<pre>
<code>ts = import_from_list(ts_list, header=TRUE) 
confounds = import_from_list(confound_list, header=FALSE)</code></pre>

<p>Next scrub the bad time points (excessive head motion) from the time series files. If you have dealt with head motion in some other way, you can skip this step.</p>

<pre>
<code>new_ts = scrub_time_series(ts, confounds, method="censor", n=100, n.nodes=116)</code></pre>

<p><br />
<strong>Connectivity Matrix</strong></p>

<p>Network construction in this package is done by creating signed, weighted networks. Allowing connections between brain regions to be positive or negative reflects the nature of functional connectivity. However, several methods are offered for doing so and which is best might depend on one&#39;s particular needs. For methods &quot;pearson&quot;, &quot;partial&quot;, &quot;covar&quot;, and &quot;parvar&quot; the estimation is done using the James-Stein based methods described in Schaefer and Strimmer (2005) and the optimal shrinkage level calculated for each subject using the method of Opgen-Rhein and Strimmer (2007). Respectively these give the regularized pearson correlation matrix, shrinkage estimated partial correlation matrix, and regularized covariance matrix. This depends on the corpcor package.The covariance option shouldn&#39;t be used directly for network analysis, but can be used with whatever you need it for. </p>


<p>After selecting your method, simply call the get_cor_matrices function.<br />
&nbsp;</p>

<pre>
<code>cormats = get_cor_matrices(new_ts, method="partial", n=100, n.nodes=116, pre.thresh=TRUE)</code></pre>

<p><br />
Several functions in rsfcNet are dependent upon or are modifications of igraph package functions. To make life easier, the mats_as_graphs will convert the correlation matrices into the igraph format.<br />
&nbsp;</p>

<pre>
<code>graphs =  mats_as_graphs(cormats)</code></pre>

<p><br />
<em><strong>Graph Theory Measures&nbsp;</strong></em><br />
<br />
The next step is selecting which graph theory measures you&#39;d like to use.&nbsp; Below is a dictionary of graph theory measures and concepts utilized in the rsfcNet package.&nbsp;</p>

<p><strong>Betweenness Centrality</strong></p>

<p>Nodes with a high betweenness centrality are considered important because they comprise pathways which control information flow. Betweenness centrality shows which nodes act as &quot;relays&quot; or &quot;hubs&quot; between nodes in a network by measuring number of times a node lies on the shortest path between other nodes. Betweenness centrality is considered a more global measure of centrality because it attempts to estimate a node&#39;s importance to network integrity. Although a popular measure in resting state fMRI research (see for example Wang, Zuo, &amp; He, 2010), its reliance on shortest paths makes it an inappropriate measure for systems like the brain which uses a parallel distributed means of information flow (where information does not deterministically follow a shortest path) (Borgatti, 2005; Telesford et al., 2011; Power et al., 2013)&nbsp;.<br />
&nbsp;</p>

<pre>
<code>For a single graph:
btwn_centr(graph) 

For multiple graphs: 
btwn_centr_mult(graphs)</code></pre>

<p><strong>Closeness Centrality</strong></p>

<p>Closeness Centrality is a measure where each node&rsquo;s importance is determined by closeness to all other nodes. Closeness is defined as the reciprocal of the sum of shortest path lengths (the distance) between <span class="math-tex">\( node_i \)</span>&nbsp;and all other nodes.</p>

<p><br />
<span class="math-tex">\(\frac{1}{∑_{j\neq i}^N d(n_i,e_{ij})}\)</span><br />
&nbsp;</p>

<p>Closeness centrality estimates how fast the flow of information would be through a given node to other nodes. Nodes with high closeness centrality may have better access to information at other nodes or more direct influence on other nodes. Although it relies on the calculation of shortest paths, it is more applicable to parallel diffusion networks than betweenness centrality. It is still apt for parallel diffusion networks because closeness centrality is measuring closeness to *any* other node, rather than treating the node as a junction between two nodes as in betweenness centrality (Borgatti, 2005).&nbsp;</p>

<p>However, see current centrality for an analagous spectral measure even more appropriate to parallel diffusion networks.</p>

<pre>
<code>For a single graph:
closeness_centr(graph)

For multiple graphs: 
closeness_centr_mult(graphs)</code></pre>

<p><strong>Clustering Coefficient (Local)&nbsp;</strong></p>

<p>The local transitivity (or clustering coeffecient) of a node is the fraction of edges a node forms with its neighbors out of the the number of edges it would take to make complete triangles. The basic method in this package is available in the local_trans function and works for weighted or unweighted networks (Barratt et al., 2004). Barrat&#39;s method is given in the formula below. For unweighted networks this simply reduces to the unweighted clustering coefficient.&nbsp;<br />
<br />
<span class="math-tex">\(C_i^w=\frac{1}{s_i(k_i-1)}\sum_{j,h}\frac{w_{ij}+w_{ih}}{2}a_{ij}a_{ih}a_{jh}\)</span></p>

<p><br />
An alternative to Barratt&#39;s method is the&nbsp;Zhang &amp; Horvath (2005) method available in the clustcoef_signed function.</p>

<p><span class="math-tex">\(C_i^w=\frac{\displaystyle\sum\nolimits_{j,q} (w_{(j,i)}w_{(i,q)}w_{(j,q)}) } {{\big({\displaystyle\sum\nolimits_{q}w_{i,q}})}^2\big)-{\displaystyle\sum\nolimits_{q}{w^2_{(i,q)}}}}.\)</span></p>

<p>&nbsp;</p>

<p>For&nbsp;<em>signed&nbsp;</em>weighted graphs one should use&nbsp;the Constantini &amp; Perugini (2014) method available in the clustcoef_signed function.</p>

<p><span class="math-tex">\(\hat{C}_{i, Z}=\frac{\sum_{j, q}\left(w_{s(j, i)} w_{s(i, q)} w_{s(j, q)}\right)}{\sum_{j \neq q}\left|w_{s(j, i)} w_{s(i, q)}\right|})</span></p>

<p>&nbsp;</p>

<pre>
<code>For a single graph (Zhang or Constantini):

clustcoef_signed(graph)

For multiple graphs (Barratt's method):

local_trans(graphs)

For multiple graphs (Zhang or Constantini):

clustcoef_signed_mult(graphs)</code></pre>

<p><br />
<strong>Clustering Coefficient (Global)</strong><br />
<br />
Simply the average of all local clustering coefficients, this tells you how closely on average connected nodes tend to be to their neighbors (Watts &amp; Strogatz, 1998)<br />
&nbsp;</p>

<pre>
<code>For a single graph:
global_clust(graph)

For multiple graphs:
global_clust_mult(graphs)

For multiple graphs (Barratt's method):
global_trans(graphs) 	</code></pre>

<p><br />
<strong>Communities ; Community Structure</strong></p>

<p>See Modules.</p>

<p><strong>Current Centrality (Circuit-Flow Closeness Centrality)</strong></p>

<p>Circuit flow closeness centrality, or current centrality, is a method of calculating the closeness more appropriate to networks where information does not travel in a necessarily serial fashion (Brandes et al., 2005). In an analogy to electricity flowing through a series of circuits, closeness is defined as the conductance between two nodes, which is the inverse of the resistance distance. Current centrality is then the average conductance between two nodes. Current centrality is a spectral measure calculated using the Moore-Penrose inversion of the graph Laplacian. Current centrality is also known as information centrality.&nbsp;</p>

<pre>
<code>For a single graph:
current_centr(graph)

For multiple graphs:
current_centr_mult(graphs)</code></pre>

<p>&nbsp;</p>

<p><strong>Diversity &nbsp;Coefficient</strong></p>

<p>A measure that characterizes the degree to which a node is connected to the entire network (high Diversity &nbsp;coefficient) or only within a module (low Diversity &nbsp;coefficient). Low Diversity &nbsp;nodes with high within-module z scores are considered provincial hubs (important for within-module communication) and high Diversity &nbsp;nodes with high within-module z-scores are considered connector hubs (important for inter-module communication).</p>

<p>The diversity coefficient is given by the following formula:</p>

<p><span class="math-tex">\(h_i = -\frac{1}{log(m)}  \sum_{u=1}^{N_M} \Bigg( \left ( \frac{s_{iu}}{s_i} \right ) \cdot log  \left ( \frac{s_{iu}}{s_i} \right ) \Bigg)\)</span><br />
<br />
Rubinov &amp; Sporns (2011) give a generalization of this to signed networks, where the diversity is calculated separately for positive and negative edge weights. These are then combined by the following formula, where&nbsp;&nbsp;<span class="math-tex">\(s_i^{-}\)</span>and <span class="math-tex">\(s_i^{+}\)</span>&nbsp;are respectively the strength of the negative and positive connections in the network. This is offered in the rsfcNet package.</p>

<p><br />
<span class="math-tex">\(h_{i}^{*} = h_{i}^{+} - \Bigg( \frac{s_i^{-}}{s_i^{+}+s_i^{-}}\Bigg) h_{i}^{-} \)</span><br />
<br />
Available as part of the module connectivity function.</p>

<pre>
<code>module_connectivity(graph, module)</code></pre>

<p><br />
<strong>Degree Centrality</strong></p>

<p>Degree centrality, also known as just degree, is simply a count of the number of non-zero connections a node has (Fornito et al., 2016b). If the nodes in a network have weighted edges, the weighted degree, also called strength, is simply the sum of a node&#39;s edge weights. Degree and Strength are considered local measures of centrality because only a given node&#39;s immediate connections are considered.<br />
<br />
For signed weighted graphs the weighted strength* measure is available:&nbsp;<br />
<br />
&nbsp;</p>

<pre>
<code>For a single graph:
degree_centr(graph)

For multiple graphs:
degree_mult(graphs)

For a weighted graph:
strength_signed(graph)

For multiple weighted graphs:
strength_multiple(graphs)</code></pre>

<p><br />
<strong>Delta Centrality (Energy)</strong></p>

<p>Delta centrality measures the change in a global property of the graph that occurs due to the deletion of a node or edge (Fornito et al., 2016a). Implemented in this package is delta energy, which tracks the change in graph energy due to each of the ith nodes being deleted. See Graph Energy for more information. Also see Laplacian Centrality and Vitality, another delta centrality measure.</p>

<pre>
<code>For a single graph:
delta_energy(graph)

For multiple graphs:
delta_energy_mult(graphs)</code></pre>

<p><br />
<strong>Neighbor centrality</strong></p>

<p>Neighbor centrality is a measure of the average degree or strength of the edges of a node&#39;s neighbors. Neighbor centrality shows which nodes are connected to well connected nodes (Barratt et al., 2004; Fornito et al., 2016). This offers an improvement over degree since a low-degree node with connections to high degree nodes may have a central role in the network. Like leverage centrality and Laplacian centrality it considers not only the immediate environment of a node but an intermediate space between the local neighborhood and global embeddedness. However, this is conceptually distinct from leverage centrality, which defines importance as being connected to nodes with fewer connections of their own.&nbsp;</p>

<pre>
<code>For a single graph:
neighbor_centr(graph)

For multiple graphs:
neighbor_centr_mult(graphs)</code></pre>

<p><br />
<strong>Eigenvector Centrality</strong></p>

<p>The eigenvector centrality is the ith entry (for the ith node) in the principal eigenvector, that is, the eigenvector belonging to the largest eigenvalue of a network (Fornito et al., 2016a). Eigenvector centrality differs conceptually from degree or strength. A node with many connections does not necessarily have a high eigenvector centrality. For example, a node may have many very weak connections that yield a large value for strength/degree. Likewise, a node with high eigenvector centrality may have a low degree but be well connected to a small number of important nodes. Eigenvector centrality is a spectral measure appropriate for networks where parallel diffusion occurs.</p>

<pre>
<code>For a single graph:
eigen_centr(graph)

For multiple graphs:
eigen_centr_mult(graphs)</code></pre>

<p><br />
<strong>Fiedler Value</strong></p>

<p>The Fiedler value is the second smallest eigenvalue of the Laplacian representation of a graph. The closer the Fiedler value is to zero the more easily the graph can be split into separate components unconnected to each other. The Fiedler value is also known as the algebraic connectivity of a graph (Mohar, 1991). Hence the Fiedler value can be used as a measure of a network&#39;s robustness to becoming disconnected.</p>

<pre>
<code>For a single graph:
fiedler_value(graph)

For multiple graphs:
fiedler_value_mult(graphs)</code></pre>

<p><br />
<strong>Graph Energy</strong></p>

<p>Graph energy was originally applied in organic chemistry to quantify the stability of molecular orbitals associated with pi-electrons (Li, Shi, &amp; Gutman, 2012). The graph energy informs about the connectivity of the graph as a whole, indicating how resilient the network might be to attack (Shatto &amp; Cetinkaya, 2017). A graph energy of zero means the nodes are not connected at all. The graph energy is calculated simply by summing the absolute values of the eigenvalues of a matrix:&nbsp;</p>

<p><br />
<span class="math-tex">\(E(G) = \sum{|\lambda_i|}\)</span><br />
&nbsp;</p>

<pre>
<code>For a single graph:
graph_energy(graph)

For multiple graphs:
graph_energy_mult(graphs)</code></pre>

<p><strong>Laplacian Centrality</strong></p>

<p>Laplacian centrality is a spectral graph theory measure and a member of the delta centrality family (where centrality is defined as the change in some graph-level measure due to the deletion of a node). Here the graph level measure of interest is the Laplacian energy of a graph, which is defined as the sum of squared eigenvalues of the graph Laplacian. This measure not only takes into account the local environment immediately around it but also the larger environment around its neighbors. It is an intermediate between metrics that assess a node&#39;s position in the whole network (such as eigenvector centrality) and the local neighborhood (such as strength).</p>

<p>The Laplacian energy can be calculated as the sum of the sums of squared degrees (weighted degree or binary) for each node and twice the sum of squared edge weights for each edge in a graph.&nbsp;</p>

<p><br />
<span class="math-tex">\(E_{L}(G)=∑_{i=1}^nd_i^2+2∑_{i&lt;j}w_{i\text{,}j}^2\)</span><br />
&nbsp;</p>

<p>The Laplacian centrality for a <span class="math-tex">\( node_i \)</span>&nbsp; is then the difference in Laplacian graph energy between the full graph and the graph where <span class="math-tex">\( node_i \)</span>&nbsp; is deleted.&nbsp;</p>

<p><br />
<span class="math-tex">\(\Delta E_{L} (G) = E_L (G) - E_L(G_{-node_i})\)</span><br />
&nbsp;</p>

<pre>
<code>For a single graph:
laplace_centr(graph)

For multiple graphs:
laplace_centr_mult(graphs)</code></pre>

<p>&nbsp;</p>

<p><strong>Leverage Centrality</strong></p>

<p>Leverage centrality defines importance as being connected to other nodes who in turn have only fewer connections.</p>

<p><span class="math-tex">\(l_i = \frac{1}{k_i} \sum_{j \in N_i} \frac{k_i - k_j}{k_i + k_j}\)</span></p>

<p>It was proposed by Joyce et al (2010) and inspired by the fact that neural connections integrate information from their connections. If a single cell (or region) synapses onto many areas that receive relatively fewer inputs, the source region will have a greater influence on the targets. Leverage centrality does not assume information flows strictly along shortest paths in a serial manner, but rather that information diffuses along parallel routes. Like neighbor centrality (its opposite) and Laplacian centrality, it can be considered as an intermediate (as opposed to global or strictly local) measure of centrality.</p>

<pre>
<code>For a single graph:
leverage_centr(graph)

For multiple graphs:
leverage_centr_mult(graphs)</code></pre>

<p><strong>Modules&nbsp;</strong></p>

<p>Modules, also known as communities, are groups of nodes that connect more (or more strongly, if considering weighted networks) to each other than to other nodes. Modules can be found by a large number of algorithms, each with its own strength and weaknesses. Module algorithms typically use clustering algorithms to find a community structure that maximizes the modularity statistic. Note that the concept of modules in the network theory sense differs from the concept in cognitive psychology of modules as functionally encapsulated regions of brain tissue. The key difference is that modules in the Fodorian sense are typically thought to be associated with a single task, while modules in the network theory sense are descriptions of how different nodes cluster together and makes no assumptions regarding *functional encapsulation*. A graph theoretic module may participate in a large number of tasks, and a single task may involve multiple modules (Colombo, 2013; Sporns &amp; Betzel, 2016).<br />
<br />
A number of methods are included in this package including the fast-greedy method (Clauset, Newman, &amp; Moore, 2004), louvain method (Blondel et al., 2008), eigenvector method (Newman, 2006), walktrap method (Gates et al., 2016; Yang et al., 2016), and re-iterative spinglass method (Traag &amp; Bruggeman, 2009; Wang et al., 2013; Zhang &amp; Moore, 2014).&nbsp;The louvain and walktrap methods are typically the better options, and run quickly as well&nbsp;&nbsp;(Gates et al., 2016; Yang et al., 2016). The fast greedy method can work well when the size of modules is not very small.&nbsp;(Yang, Algesheimer, &amp; Tessone, 2016)</p>

<pre>
<code>For a single graph: 
get_modules(graph, method="louvain")

Applied to a list of graphs:
lapply(graphs, function(i) get_modules(i, method="louvain"))
</code></pre>

<p><strong>Participation Coefficient</strong></p>

<p>A measure that characterizes the degree to which a node is connected to the entire network (high participation coefficient) or only within a module (low participation coefficient) (Guimera &amp; Nunes Amaral, 2005). Low participation nodes with high within-module z scores are considered provincial hubs (important for within-module communication) and high participation nodes with high within-module z-scores are considered connector hubs (important for inter-module communication).</p>

<p>The participation coefficient is defined by the following formula:<br />
<br />
<span class="math-tex">\(p_i = 1 - \sum_{s=1}^{N_M} \left ( \frac{s_{iu}}{s_i} \right )^2\)</span></p>

<p>One can define a weighted participation for signed networks analogous to the diversity coefficient, which is offered in the rsfcNet package:</p>

<p><span class="math-tex">\(p_{i}^{*} = p_{i}^{+} - \Bigg( \frac{s_i^{-}}{s_i^{+}+s_i^{-}} \Bigg) p_{i}^{-} \)</span></p>

<pre>
<code>module_connectivity(graph, module)</code></pre>

<p><br />
<strong>Strength&nbsp;</strong></p>

<p>See Degree Centrality.<br />
<br />
<strong>Transitivity (Local)</strong><br />
<br />
See clustering coefficient (local).<br />
<br />
<strong>Transitivity (Global)</strong><br />
<br />
See clustering coefficient (global)&nbsp;</p>

<p><strong>Vitality&nbsp;</strong></p>

<p>Vitality, also known as closeness vitality, is a delta measure of centrality. The closeness vitality of a node is the change in total distance between all other nodes in a graph when a node is deleted (Brandes &amp; Erlebach, 2005). A node with high closeness vitality creates greater distance between nodes in the network when it is deleted, implying it has a privileged place in the network and is vital to global communication.&nbsp;</p>

<pre>
<code>For a single graph:
vitality(graph)

For multiple graphs:
vitality(graphs)</code></pre>


<strong><p>Bibliography</p></strong>

<p>Barrat, A., Barth&eacute;lemy, M., Pastor-Satorras, R., &amp; Vespignani, A. (2004). The architecture of complex weighted networks. <em>Proceedings of the National Academy of Sciences of the United States of America</em>, <em>101</em>(11), 3747&ndash;3752. <a href="https://doi.org/10.1073/pnas.0400087101">https://doi.org/10.1073/pnas.0400087101</a></p>

<p>Blondel, V. D., Guillaume, J.-L., Lambiotte, R., &amp; Lefebvre, E. (2008). Fast unfolding of communities in large networks. <em>Journal of Statistical Mechanics: Theory and Experiment</em>, <em>2008</em>(10), P10008. <a href="https://doi.org/10.1088/1742-5468/2008/10/P10008">https://doi.org/10.1088/1742-5468/2008/10/P10008</a></p>

<p>Borgatti, S. P. (2005). Centrality and network flow. <em>Social Networks</em>, <em>27</em>(1), 55&ndash;71. <a href="https://doi.org/10.1016/j.socnet.2004.11.008">https://doi.org/10.1016/j.socnet.2004.11.008</a></p>

<p>Brandes, U., &amp; Erlebach, T. (Eds.). (2005). <em>Network Analysis</em> (Vol. 3418). Berlin, Heidelberg: Springer Berlin Heidelberg. <a href="https://doi.org/10.1007/b106453">https://doi.org/10.1007/b106453</a></p>

<p>Brandes, U., &amp; Fleischer, D. (2005). Centrality Measures Based on Current Flow. In V. Diekert &amp; B. Durand (Eds.), <em>STACS 2005</em> (pp. 533&ndash;544). Springer Berlin Heidelberg.</p>

<p>Clauset, A., Newman, M. E. J., &amp; Moore, C. (2004). Finding community structure in very large networks. <em>Physical Review E</em>, <em>70</em>(6). <a href="https://doi.org/10.1103/PhysRevE.70.066111">https://doi.org/10.1103/PhysRevE.70.066111</a></p>

<p>Colombo, M. (2013). Moving Forward (and Beyond) the Modularity Debate: A Network Perspective. <em>Philosophy of Science</em>, <em>80</em>(3), 356&ndash;377. <a href="https://doi.org/10.1086/670331">https://doi.org/10.1086/670331</a></p>

<p>Costantini, G., &amp; Perugini, M. (2014). Generalization of Clustering Coefficients to Signed Correlation Networks. <em>PLoS ONE</em>, <em>9</em>(2), e88669. <a href="https://doi.org/10.1371/journal.pone.0088669">https://doi.org/10.1371/journal.pone.0088669</a></p>

<p>Csardi, G., &amp; Nepusz, T. (2006). The igraph software package for complex network research. <em>InterJournal</em>, <em>Complex Systems</em>, 1695.</p>

<p>Fornito, A., Zalesky, A., &amp; Bullmore, E. (2016a). Centrality and Hubs. In <em>Fundamentals of Brain Network Analysis</em> (pp. 137&ndash;161). Elsevier. <a href="https://doi.org/10.1016/B978-0-12-407908-3.00005-4">https://doi.org/10.1016/B978-0-12-407908-3.00005-4</a></p>

<p>Fornito, A., Zalesky, A., &amp; Bullmore, E. (2016b). Node Degree and Strength. In <em>Fundamentals of Brain Network Analysis</em> (pp. 115&ndash;136). Elsevier. <a href="https://doi.org/10.1016/B978-0-12-407908-3.00004-2">https://doi.org/10.1016/B978-0-12-407908-3.00004-2</a></p>

<p>Gates, K. M., Henry, T., Steinley, D., &amp; Fair, D. A. (2016). A Monte Carlo Evaluation of Weighted Community Detection Algorithms. <em>Frontiers in Neuroinformatics</em>, <em>10</em>. <a href="https://doi.org/10.3389/fninf.2016.00045">https://doi.org/10.3389/fninf.2016.00045</a></p>

<p>Guimer&agrave;, R., &amp; Nunes Amaral, L. A. (2005). Functional cartography of complex metabolic networks. <em>Nature</em>, <em>433</em>(7028), 895&ndash;900. <a href="https://doi.org/10.1038/nature03288">https://doi.org/10.1038/nature03288</a></p>

<p>Joyce, K. E., Laurienti, P. J., Burdette, J. H., &amp; Hayasaka, S. (2010). A New Measure of Centrality for Brain Networks. <em>PLoS ONE</em>, <em>5</em>(8), e12200. <a href="https://doi.org/10.1371/journal.pone.0012200">https://doi.org/10.1371/journal.pone.0012200</a></p>

<p>Kr&auml;mer, N., Sch&auml;fer, J., &amp; Boulesteix, A.-L. (2009). Regularized estimation of large-scale gene association networks using graphical Gaussian models. <em>BMC Bioinformatics</em>, <em>10</em>(1), 384. <a href="https://doi.org/10.1186/1471-2105-10-384">https://doi.org/10.1186/1471-2105-10-384</a></p>

<p>Li, X., Shi, Y., &amp; Gutman, I. (2012). <em>Graph Energy</em>. New York, NY: Springer New York. <a href="https://doi.org/10.1007/978-1-4614-4220-2">https://doi.org/10.1007/978-1-4614-4220-2</a></p>

<p>Lohmann, G., Margulies, D. S., Horstmann, A., Pleger, B., Lepsien, J., Goldhahn, D., &hellip; Turner, R. (2010). Eigenvector Centrality Mapping for Analyzing Connectivity Patterns in fMRI Data of the Human Brain. <em>PLoS ONE</em>, <em>5</em>(4), e10232. <a href="https://doi.org/10.1371/journal.pone.0010232">https://doi.org/10.1371/journal.pone.0010232</a></p>

<p>Mohar, B. (1991). Eigenvalues, diameter, and mean distance in graphs. <em>Graphs and Combinatorics</em>, <em>7</em>(1), 53&ndash;64. <a href="https://doi.org/10.1007/BF01789463">https://doi.org/10.1007/BF01789463</a></p>

<p>Newman, M. E. J. (2006). Finding community structure in networks using the eigenvectors of matrices. <em>Physical Review E</em>, <em>74</em>(3). <a href="https://doi.org/10.1103/PhysRevE.74.036104">https://doi.org/10.1103/PhysRevE.74.036104</a></p>

<p>Opgen-Rhein, R., &amp; Strimmer, K. (2007). Accurate Ranking of Differentially Expressed Genes by a Distribution-Free Shrinkage Approach. <em>Statistical Applications in Genetics and Molecular Biology</em>, <em>6</em>(1). <a href="https://doi.org/10.2202/1544-6115.1252">https://doi.org/10.2202/1544-6115.1252</a></p>

<p>Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., &amp; Petersen, S. E. (2012). Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. <em>NeuroImage</em>, <em>59</em>(3), 2142&ndash;2154. <a href="https://doi.org/10.1016/j.neuroimage.2011.10.018">https://doi.org/10.1016/j.neuroimage.2011.10.018</a></p>

<p>Power, J. D., Schlaggar, B. L., Lessov-Schlaggar, C. N., &amp; Petersen, S. E. (2013). Evidence for Hubs in Human Functional Brain Networks. <em>Neuron</em>, <em>79</em>(4), 798&ndash;813. <a href="https://doi.org/10.1016/j.neuron.2013.07.035">https://doi.org/10.1016/j.neuron.2013.07.035</a></p>

<p>Qi, X., Fuller, E., Wu, Q., Wu, Y., &amp; Zhang, C.-Q. (2012). Laplacian centrality: A new centrality measure for weighted networks. <em>Intelligent Knowledge-Based Models and Methodologies for Complex Information Systems</em>, <em>194</em>, 240&ndash;253. <a href="https://doi.org/10.1016/j.ins.2011.12.027">https://doi.org/10.1016/j.ins.2011.12.027</a></p>

<p>Rubinov, M., &amp; Sporns, O. (2011). Weight-conserving characterization of complex functional brain networks. <em>Neuroimage</em>, <em>56</em>(4), 2068&ndash;2079. <a href="https://doi.org/10.1016/j.neuroimage.2011.03.069">https://doi.org/10.1016/j.neuroimage.2011.03.069</a></p>

<p>Sch&auml;fer, J., &amp; Strimmer, K. (2005). A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics. <em>Statistical Applications in Genetics and Molecular Biology</em>, <em>4</em>(1). <a href="https://doi.org/10.2202/1544-6115.1175">https://doi.org/10.2202/1544-6115.1175</a></p>

<p>Shatto, T. A., &amp; Cetinkaya, E. K. (2017). Variations in graph energy: A measure for network resilience. In <em>2017 9th International Workshop on Resilient Networks Design and Modeling (RNDM)</em> (pp. 1&ndash;7). Alghero, Italy: IEEE. <a href="https://doi.org/10.1109/RNDM.2017.8093019">https://doi.org/10.1109/RNDM.2017.8093019</a></p>

<p>Sporns, O., &amp; Betzel, R. F. (2016). Modular Brain Networks. <em>Annual Review of Psychology</em>, <em>67</em>, 613&ndash;640. <a href="https://doi.org/10.1146/annurev-psych-122414-033634">https://doi.org/10.1146/annurev-psych-122414-033634</a></p>

<p>Telesford, Q. K., Simpson, S. L., Burdette, J. H., Hayasaka, S., &amp; Laurienti, P. J. (2011). The Brain as a Complex System: Using Network Science as a Tool for Understanding the Brain. <em>Brain Connectivity</em>, <em>1</em>(4), 295&ndash;308. <a href="https://doi.org/10.1089/brain.2011.0055">https://doi.org/10.1089/brain.2011.0055</a></p>

<p>Traag, V. A., &amp; Bruggeman, J. (2009). Community detection in networks with positive and negative links. <em>Physical Review E</em>, <em>80</em>(3). <a href="https://doi.org/10.1103/PhysRevE.80.036115">https://doi.org/10.1103/PhysRevE.80.036115</a></p>

<p>Wang, J., Zuo, X., &amp; He, Y. (2010). Graph-based network analysis of resting-state functional MRI. <em>Frontiers in Systems Neuroscience</em>. <a href="https://doi.org/10.3389/fnsys.2010.00016">https://doi.org/10.3389/fnsys.2010.00016</a></p>

<p>Wang, Z., Hu, Y., Xiao, W., &amp; Ge, B. (2013). Overlapping community detection using a generative model for networks. <em>Physica A: Statistical Mechanics and Its Applications</em>, <em>392</em>(20), 5218&ndash;5230. <a href="https://doi.org/10.1016/j.physa.2013.06.038">https://doi.org/10.1016/j.physa.2013.06.038</a></p>

<p>Watts, D. J., &amp; Strogatz, S. H. (1998). Collective dynamics of &lsquo;small-world&rsquo; networks. <em>Nature</em>, <em>393</em>(6684), 440&ndash;442. <a href="https://doi.org/10.1038/30918">https://doi.org/10.1038/30918</a></p>

<p>Yang, Z., Algesheimer, R., &amp; Tessone, C. J. (2016). A Comparative Analysis of Community Detection Algorithms on Artificial Networks. <em>Scientific Reports</em>, <em>6</em>(1). <a href="https://doi.org/10.1038/srep30750">https://doi.org/10.1038/srep30750</a></p>

<p>Zhang, B, &amp; Horvath, S. (2005). A General Framework for Weighted Gene Co-Expression Network Analysis. <em>Statistical Applications in Genetics and Molecular Biology</em>, <em>4</em>(1). <a href="https://doi.org/10.2202/1544-6115.1128">https://doi.org/10.2202/1544-6115.1128</a></p>

<p>Zhang, P., &amp; Moore, C. (2014). Scalable detection of statistically significant communities and hierarchies, using message passing for modularity. <em>Proceedings of the National Academy of Sciences</em>, <em>111</em>(51), 18144&ndash;18149. <a href="https://doi.org/10.1073/pnas.1409770111">https://doi.org/10.1073/pnas.1409770111</a></p>

<p>Zhu, Y., &amp; Cribben, I. (2018). Sparse Graphical Models for Functional Connectivity Networks: Best Methods and the Autocorrelation Issue. <em>Brain Connectivity</em>, <em>8</em>(3), 139&ndash;165. <a href="https://doi.org/10.1089/brain.2017.0511">https://doi.org/10.1089/brain.2017.0511</a></p>

<p>Zou, H. (2006). The Adaptive Lasso and Its Oracle Properties. <em>Journal of the American Statistical Association</em>, <em>101</em>(476), 1418&ndash;1429. <a href="https://doi.org/10.1198/016214506000000735">https://doi.org/10.1198/016214506000000735</a></p>

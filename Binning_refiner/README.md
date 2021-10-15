# Binning refiner

This project is based on the metabolite connectivity score idea from this [paper](https://academic.oup.com/bioinformatics/article/32/6/867/1744247). In this paper they used the connectivity of the metabolic network to bin unbinned contigs.
Unfortunately, this did not work as the metabolite connectivity score did not provide enough evidence to be used as a basis for binning refinement (additional binning of highly complete bins) . We tried to add other features to help select potential candidate contigs (i.e., contigs to assign to bins), unfortunately since the MCS did not produce a strong signal, we decided to abandon this project.

If you want to have access to the graph database please contact me.

input files:

- annotation_CDS_RNA_hmms_checkm.gff - prokka gff for masking certain regions
- assembly.fa
- bins.tsv -  first column should be contig ID and second the bin id
- consensus_annotation.tsv - annotation from Mantis using uniprot_ec and uniprot_rhea as reference - check https://github.com/PedroMTQ/mantis and https://github.com/PedroMTQ/refdb_generator
- contig_depth.txt - first column should be contig ID and second the depth
- taxonomy_lineage.tsv - `contig	superkingdom	phylum	class	order	family	genus	species`




### Workflow

Our method can be summarized in five major steps: (i) network generation, (ii) creation of reference HMMs, (iii) pre-processing, (iv) metabolic information extraction, and (v) contigs assignment. Steps one and two are used for metabolic information retrieval and functional annotation, respectively; these have been previously executed and the associated files are available for download. Steps three to five are executed during the binning refinement.

#### Network generation

To create the global metabolic network we started by downloading and parsing the Rhea and ChEBI databases, this data was then structured into a metabolic network (i.e., compounds, reactions and proteins) and saved into a Neo4j database. We then compiled a list of all KEGG KOs and collected information on these KOs from KEGG, linking these to their respective proteins, reactions and compounds; this data was then added to the previously compiled Neo4j database by cross-linking the information from Rhea and KEGG, thus creating a highly comprehensive metabolic network. In total, the Neo4j metabolic network database contained a total of 299742 nodes, 6349 protein nodes, 20186 reaction nodes, and 25172 compound nodes (the remaining nodes correspond to information linked to the previous entities, e.g., database identifiers). This Neo4j database is available at XXX.


#### Functional annotation reference HMMs

In order to functionally annotate the assembly, and thus understand its metabolic potential, we used Mantis. In order to be able to cross-link our Neo4j metabolic network we used Mantis with the KOfam database (one of its default reference datasets), and additionally compiled two reference datasets. Specifically, we compiled two HMMs, one clustered via Rhea IDs and another via EC IDs. To compile these HMMs we download all SwissProt (CITE) protein sequences; for the EC-centric HMMs we downloaded SwissProt's protein sequences metadata and generated one fasta per EC; For the Rhea-centric HMMs we downloaded Rhea's metadata (i.e., rhea2uniprot and rhea2xrefs) and generated one fasta per Rhea ID and also a Mantis-compatible metadata file to allow for cross-linking of Rhea IDs with other IDs). In both cases, we defined a minimum of 10 sequences per fasta. We then applied the following workflow to generate both HMMs: (i) create a multiple sequence alignment (using MUSCLE (CITE) for fastas with 500 or less sequences, and CLUSTALO (CITE) for fastas with more than 500 sequences), (ii) generate and press HMMs profiles with HMMER (iii) merge HMMs into one global HMM. The code for this process is available at https://github.com/PedroMTQ/hmm\_updater

#### Pre-processing

The pre-processing is responsible for collecting data on the assembly (i.e., contig kmer frequency, depth and taxonomy), its respective bins, and the functional annotations for each of the assembly's contigs CDS.
Pre-processing contains several steps: (i) identification of contig features (i.e. CDS and low complexity features, e.g., rRNA or CRISPR regions); (ii) parsing of Mantis-derived functional annotations; (iii) assembly parsing where valid contigs (i.e., a minimum bp length of 300 if it contains a functional annotation or 1000 if it does not) are identified and extracted. During step (iv), low complexity features are masked so that during during step (v) the kmer frequency counting (for kmer sizes of 2,3, and 4) does not include these regions. This step results in a matrix containing kmer frequency per contig, which is then scaled in step (vi) according to a previously published methodology by Hickl et al. (BINNY). In step (vii), contig depth and taxonomy is added to this matrix and in step (viii) this matrix is stored in a SQLite database for later retrieval. In step (ix) binned and unbinned contigs are extracted from the binny-derived binning output, and maximum intra-bin contig distance is calculated by calculating the cosine similarity (1-cosine distance) between each contig in the bin to all other contigs in the same bin (this is applied to every bin).
We also used paired-end reads to assign unassigned contigs to bins.


#### Metabolic information extraction

During metabolic information extraction, functional information from the assembly is linked with the Neo4j database.
During this step we start by calculating the weight of each compound within the Neo4j network (N=25172) by applying a technique similar to inverse document frequency (IDF), a technique commonly used in text mining (CITATION). This was done by checking how many connections ${CP}$ each compound node has to all reaction nodes (e.g., if the compound water appears in 3 reactions in the whole network, then ${CP}$ for water is equal to 3) and dividing the amount of reactions ${AR}$ in the network by ${CP}$, thus ${IDF=\frac{AR}{CP}}$. This method gives a higher weight to compounds that appear less often and less weight to compounds that appear more often (i.e., more likely to be cofactors). This weight was then scaled by XXXXXXXXXXXX so that cofactors or more common compounds still remain representative during the metabolite connectivity score (${MCS}$) calculation. Information on the binned and unbinned reactions is then collected from the Neo4j database (i.e., reactants/products and reaction reversibility). Each unbinned contig will then have a set of reactions associated with it, as well as a set of reactants and products per reaction; each bin will have a set of of metabolites that the bin's organism should produce or consume (as derived from the functional annotations).

\subsubsection{Unbinned contigs assignments}

%https://academic.oup.com/bioinformatics/article/32/6/867/1744247
In this last step, we measure the ${MCS}$ of each bin per reaction, where ${MCS}$ is defined as:
\[
\frac{R1+...+Rn}{NR} + \frac{P1+...+Pn}{NP}
\]
Where R stands for reactant, NR the number of reactants, P for products and NP the number of products.
For example, if a contig contains a CDS functionally annotated with the EC 3.1.3.10, then by extent this contig is involved in the reaction: D-Glucose 1-phosphate + H2O <=> D-Glucose + Phosphate (two products and two reactants); We then check the ${MCS}$ per bin, where, e.g., if the candidate bin contains (within its set of produced/consumed metabolites) the metabolites Phosphate and H2O, then ${MCS=\frac{0+1}{2}+\frac{0+1}{2}}=1$. Since a contig may contain more than one reaction, we calculate the average ${MCS}$ score per bin/contig (but rejecting null ${MCS}$). It is important to note that during ${MCS}$ calculation we only assume a bin is a potential candidate if it does not contain the functional annotation being currently assessed.
We then also calculate the kmer frequency and depth distance score (${KDD}$, 0-1), and taxa distance (${TD}$,0-1) between the candidate bin and the current contig, such that the final score is equal to ${MCS*KDD*TD}$.
The kmer frequency and depth distance score equals to the average cosine similarity of each contig in the candidate bin to the current contig (using the vectors previously stored in SQLite); the taxa distance equals to the average lineage distance between each contig in the candidate bin to the current contig.

Lineage distance is calculated by adding X distance between each mismatched taxonomic rank, X being: (i) from species to genus, 0.5; (ii) genus to family, 1; (iii) family to order, 1.5; (iv) order to class, 2; (v) class to phylum, 2.5; and (vi) phylum to superkingdom, 3. For example, the distance between the \textit{Bacillus subtilis} and  \textit{Bacillus infernus} is 0.5 (different species only), but the distance between \textit{Bacillus subtilis} and \textit{Staphylococcus aureus} is 0.5 (different species) plus 1 (different genus). Since the maximum possible distance is 10.5 (0.5+1+1.5+2+2.5+3) then the taxa distance score can be calculated by dividing the taxa distance by the maximum distance and then subtracting this value to 1 (since for MCS and KDD, a higher score implies that the bin is a better candidate).  

After assigning the current contig to the best candidate bin, each bin's set of metabolites is complemented with the metabolites associated with the newly assigned contigs.

This cycle of assigning contigs is repeated until no more contigs can be assigned.


In addition to MCS, in this step we also use the kegg modules completeness (selecting the most complete kegg module pathway and checking how much more complete it gets by adding a contig)

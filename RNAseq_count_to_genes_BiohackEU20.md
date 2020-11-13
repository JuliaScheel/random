

---

**RNA-seq counts to genes for follow up in Minerva**

GEO: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz"
The Galaxy history can be found [here](https://usegalaxy.eu/u/jscheel1/h/gse147507-to-de-limma)

tags:
  - limma-voom
  - human
  - Covid19
  
questions:
  - "What are the differentially expressed genes in mock infected epithelial cells vs Covid19 infected epithelial cells?"
  - "How to analyze RNA count data using limma-voom?"
  
objectives:
  - "Analysis of RNA-seq count data using limma-voom"
  - "Visualisation and interactive exploration of count data"
  - "Identification of differentially expressed genes"
  
key_points:
  - "The limma-voom tool can be used to perform differential expression and output useful plots"
Based on: [RNA counts to genes tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.html)

---


# Introduction
The Covid dataset was integrated as part of BioHackEU20 as part of a workflow to integrate Covid19 data to a gene expression that can be visualized in Minerva and Wikpathways.

**Human Covid19 data set**

The data for this tutorial comes from a Nature Cell Biology paper, [SARS-CoV-2 launches a unique transcriptional signature from in vitro, ex vivo, and in vivo systems](https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1)), Blanco-Melo et al. 2020. The processed RNA-seq data (counts) can be downloaded from Gene Expression Omnibus database (GEO) under accession number [GSE147507](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507).

This study examined the expression profiles of  the host response to SARS-CoV-2 as it compares to other respiratory infections. Cell models of SARS-CoV-2 infection and ranscriptional profiling of a COVID-19 lung biopsy are available.

78 groups are present, with one for each combination of cell type and infection status. In this tutorial we will use the GEO counts file as a starting point for our analysis. Alternatively, you could create a count matrix from the raw sequence reads, as demonstrated in the [RNA-seq reads to counts tutorial]({% link topics/transcriptomics/tutorials/rna-seq-reads-to-counts/tutorial.md %}).

We will use **limma-voom** for identifying differentially expressed genes here. Other popular alternatives are edgeR and DESeq2. Limma-voom has been shown to be perform well in terms of precision, accuracy and sensitivity ([Costa-Silva, Domingues and Lopes 2017](https://www.ncbi.nlm.nih.gov/pubmed/29267363)) and, due to its speed, it's particularly recommended for large-scale datasets with 100s of samples ([Chen, Lun, Smyth 2016](https://f1000research.com/articles/5-1438/v2)).

This is a Galaxy tutorial based on material from the [RNA counts to genes tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.html)


# Preparing the inputs

We will use two files for this analysis:

 * **Count matrix** (genes in rows, samples in columns)
 * **Sample information** file (sample id, group)

## Import data

> ### Hands-on: Data upload
>
> 1. Create a new history for this RNA-seq exercise e.g. `RNA-seq count to genes with limma`
>
>    * Click the **+** icon at the top of the history panel
>
> 2. Import the human Covid19 counts table and the associated sample information file.
>
>     To import the files:
>     - From [GEO](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz)
>     - You can paste the link below into the **Paste/Fetch** box:
>
>     ```
>     https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts_Human.tsv.gz
>     ```
>
> 2. Rename the counts dataset as `seqdata` and the sample information dataset as `sampleinfo` using the **pencil** icon.
> 3. Check that the datatype is `tabular`.
>    If the datatype is not `tabular`, please change the file type to `tabular`.
{: .hands_on}


Let’s take a look at the data. The `seqdata` file contains information about genes (one gene per row), the first column has the Entrez gene id, the second has the gene length and the remaining columns contain information about the number of reads aligning to the gene in each experimental sample. There are two replicates for each cell type and time point (detailed sample info can be foun in GEO Accession viewer under "Overall design"). The first few rows and columns of the seqdata file are shown below.

![seqdata file](https://github.com/JuliaScheel/random/blob/main/images/countdata.PNG "Count file (before formatting)")

The `sampleinfo` file contains basic information about the samples that we will need for the analysis. See below.

![sampleinfo file](https://github.com/JuliaScheel/random/blob/main/images/sampleinfo.PNG "Sample information file (before formatting)")

## Format the data

> ### Hands-on: Format the counts data
> Rename file as `countdata` using the **pencil** icon. 

Next, let's create a new file, `factordata`, that contains the groups information that we need for the limma-voom tool. We'll combine the cell type and infection status to make 78 groups e.g. we'll combine the CellType `NHBE` with the Status `mock`. We'll use the **Merge Columns** tool to combine the cell type and mouse status columns in the sample information file, making a column with the 6 group names.

> ### Hands-on: Format the sample information file
>
> 1. **Merge Columns together** with the following parameters:
>      -  *"Select data"*: `sampleinfo`
>      -  *"Merge column"*: `Column: 3`
>      -  *"with column"*: `Column: 4`
> 2. **Cut columns from a table (cut)**  with the following parameters:
>      -  *"File to cut"*: output of **Merge Columns** 
>      -  *"Operation"*: `Keep`
>      -  *"List of fields"*: Select `Column:1` and `Column:5`
> 3. Rename file as `factordata` using the **pencil** icon. The file should look like below.
>    ![factordata file](https://github.com/JuliaScheel/random/blob/main/images/factordata.PNG "Sample information file (after formatting)")
{: .hands_on}

# Differential expression with limma-voom

## Filtering to remove lowly expressed genes

It is recommended to filter for lowly expressed genes when running the limma-voom tool. Genes with very low counts across all samples provide little evidence for differential expression and they interfere with some of the statistical approximations that are used later in the pipeline. They also add to the multiple testing burden when estimating false discovery rates, reducing power to detect differentially expressed genes. These genes should be filtered out prior to further analysis.

There are a few ways to filter out lowly expressed genes. When there are biological replicates in each group, in this case we have a sample size of 2 in each group, we favour filtering on a minimum counts-per-million (CPM) threshold present in at least 2 samples. Two represents the smallest sample size for each group in our experiment. In this dataset, we choose to retain genes if they are expressed at a CPM above 0.5 in at least two samples. The CPM threshold selected can be compared to the raw count with the CpmPlots (see below).

> ### More details on filtering
>
> The limma tool uses the `cpm` function from the edgeR package [Robinson, McCarthy, and Smyth 2010](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2796818/) to generate the CPM values which can then be filtered. Note that by converting to CPMs we are normalizing for the different sequencing depths for each sample. A CPM of 0.5 is used as it corresponds to a count of 10-15 for the library sizes in this data set. If the count is any smaller, it is considered to be very low, indicating that the associated gene is not expressed in that sample. A requirement for expression in two or more libraries is used as each group contains two replicates. This ensures that a gene will be retained if it is only expressed in one group. Smaller CPM thresholds are usually appropriate for larger libraries. As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5. You should filter with CPMs rather than filtering on the counts directly, as the latter does not account for differences in library sizes between samples.
{: .details}

## Normalization for composition bias

In an RNA-seq analysis, the counts are normalized for different sequencing depths between samples. Normalizing to eliminate composition biases between samples is also typically performed. Composition biases can occur, for example, if there are a few highly expressed genes dominating in some samples, leading to less reads from other genes. By default, TMM normalization [(Robinson and Oshlack 2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2864565/) is performed by the limma tool using the edgeR `calcNormFactors` function (this can be changed under **Advanced Options**). 

## Specify Contrast(s) of interest

Since we are interested in differences between groups, we need to specify which comparisons we want to test. For example, if we are interested in knowing which genes are differentially expressed between the mock and Covid19 infected group in the NHBE cells we specify `NHBEmock-NHBECovid` for the *Contrast of Interest*. Note that the group names in the contrast must exactly match the names of the groups in the `factordata` file. More than one contrast can be specified using the `Insert Contrast` button.

> ### Hands-on: Differential expression with limma-voom
>
> 1. **limma** with the following parameters:
>      -  *"Differential Expression Method"*: `limma-voom`
>      -  *"Count Files or Matrix?*": `Single Count Matrix`
>          -  *"Count Matrix"*: Select `countdata`
>      -  *"Input factor information from file?"*: `Yes`
>          -  *"Factor File"*: Select `factordata`
>      -  *"Use Gene Annotations?"*: `No`
>      -  *"Contrast of Interest"*: `NHBEmock-NHBECovid`
>      -  *"Filter lowly expressed genes?"*: `Yes`
>          -  *"Filter on CPM or Count values?"*: `CPM`
>          -  *"Minimum CPM"*: `0.5`
>          -  *"Minimum Samples"*: `2`
> 2. Inspect the `Report` produced by clicking on the **eye** icon
{: .hands_on}

At this point you would usually spend some time on quality control. Double check that your samples are labelled correctly, consult the multidimensional scaling plots, volcano plots, and voom variance plots created in the 'Report' according to this [tutorial](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-counts-to-genes/tutorial.html). Every step within the Galaxy history can be re-run with corrected parameters.



## Formatting results to be useable in Minerva
We require a specific format for our results to be mapped onto Minerva: Gene ID, logFC, p, adj_p
Click on the 'limma-voom_NHBEmock-NHBECovid.tsv' to download the results and upload them back into your history for further use. 

> ### Hands-on: Format the differentially expressed genes file
>
> 1. **Cut Columns from a table (cut)** with the following parameters:
>      -  *"File to cut"*: output of **limma-voom_NHBEmock-NHBECovid.tsv** 
>      -  *"Operation"*: `Keep`
>      -  *"List of fields"*: Select `Column:1`, 'Column:2, `Column:5`, 'Column:6', 'Column:4', and 'Column:7' 

> The file should look like below.
>    ![factordata file](https://github.com/JuliaScheel/random/blob/main/images/DE_sorted.PNG)
{: .hands_on}


# Conclusion

We have converted a count file into differentially expressed genes with limma-voom, which can then be visualized in Minerva.


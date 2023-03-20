## What is Cellxgene?

Cellxgene is an application developed by Celsius Tx to exclusively find, download, explore, analyze, annotate, and publish single-cell sequencing datasets.

## When should CellxGene be used?

1. **To examine categorical metadata:** Categorical metadata (such as tissue of origin or cell type) can be visualized and examined in a number of ways within cellxgene such as coloring embedding plots (i.e. color UMAP by cell type), looking at cell counts, making selections of cells or viewing the interaction between different categorical metadata fields.
2. **To find cells where a gene is expressed:** Numerical metadata about gene expression features or the number of genes can be examined on the embedding plot and be used to filter and select cells.
3. **To select and subset cells:** Cells in the embedding plot can be selected based on the gene expression cutoffs, and categorical metadata attributes.
4. **To compare the expression of multiple genes:** Cellxgene allows user to compare the expression of multiple genes via bivariate plots.
5. **To use gene sets to learn about cell population functional characteristics:** cellxgene allows user to examine groups of genes via the gene sets feature
6. **To find Marker Genes:** Cellxgene allows user to find marker genes between selected cell populations

## Cellxgene interface

The cellxgene interface is divided into 3 sections -

- Left-hand sidebar contains categorical metadata and fields curated for each dataset ingested in OmixAtlas.
- Center panel contains an embedding plot where each dot represents a cell. Cellxgene also allows users to choose different embeddings based on their needs.
- Right-hand sidebar contains numerical metadata (QC metadata) and information about gene and gene sets. The gene plots on this sidebar automatically appears for each curated dataset opened through Polly OmixAtlas.

![Homepage](../img/Cellxgene/1.png) <center>**Figure 1.** Cellxgene interface</center>

![Embed_choice](../img/Cellxgene/2.png) <center>**Figure 2.** Embedding plot options</center>

## Examining categorical metadata

Categorical data can be examined in multiple ways -

1. Colour embedding plots - The plot can be colored by different metadata. The plot can be colored by clicking on the drop icon next to metadata of interest. Color code is displayed next to the cell count of the chosen metadata.
2. Cell counts - Each metadata category will have the cell counts along with the color code.

![Example](../img/Cellxgene/3.png) <center>**Figure 3.** Example 1 - Embedding colored by cell type</center>


![Example](../img/Cellxgene/4.png) <center>**Figure 4.** Example 2 - Embedding colored by disease</center>


3. Making the selection of cells - Upon hovering over the metadata on left side, user can see the type of cells highlighted in the embedding plot.

![Highlight](../img/Cellxgene/5.png) <center>**Figure 5.** Highlighting via metadata</center>

Users can click on the 'display categories' icon to display the labels over the cluster centroid.

![display categories](../img/Cellxgene/6.png) <center>**Figure 6.** Display categories</center>

Users can select/deselect the cells by choosing the checkbox next to the value in the categorical metadat field. The unselected cells will have a smaller point size on the embedding plot.

![Highlight](../img/Cellxgene/7.png) <center>**Figure 7.** Select/deselect metadata</center>

Users can use the checkmark on the parent metadata to delect all the cells and then select the cells of interest, this makes it easier to highlight the cells that pertain to a specific tissue.

![Highlight](../img/Cellxgene/8.png) <center>**Figure 8.** Select/deselect metadata</center>

4. Viewing the relationship between different categorical fields - After users color by a particular metadata category, for example 'Cell type', users can see the distribution of the cell type in any other category by expanding that categorical metadata field, for example 'tissue'. In the plot below, cell types belonging to a particular tissue are clearly shown with the colors in the bar.

![Highlight](../img/Cellxgene/9.png) <center>**Figure 9.** Relationship between different categorical fields</center>

## Finding cells where gene is expressed

Numerical metadata on the embedding plot can be examined and used to filter and select cells. The numerical data is present on the right hand sidebar and users can click on the droplet icon to color the plot by qc metrics (for example, n\_genes\_by\_counts - number of genes that have been detected in the cell).

![Highlight](../img/Cellxgene/10.png) <center>**Figure 10.** Numerical metadata</center>

To deal with outliers, users can clip the data by clicking on 'clip' icon on top of the toolbar. For example, here we have set the values to 0 to 99 percentile and we can observe the change in the scale, graph and the embedding plot.

![Highlight](../img/Cellxgene/11.png) <center>**Figure 11.** Clipping tool</center>

To find the gene expression pattern for a particular gene, type the gene name in the search bar.

![Highlight](../img/Cellxgene/12.png) <center>**Figure 12.** Search bar for gene expression pattern</center>

For example, if we search for STAT3, a graph will appear on the top right side that contains

- Gene name
- A histogram depicting the gene expression
- Remove icon
- X-plot/Y-plot for bivariate plots
- Droplet icon to color the cells based on the gene expression

![Highlight](../img/Cellxgene/13.png) <center>**Figure 13.** Colored embedding plot via gene expression</center>

Users can use the clipping tool to remove the outliers to understand the gene expression better.

![Highlight](../img/Cellxgene/14.png) <center>**Figure 14.** Using clippng tool to study gene expression</center>

## Selecting and subsetting cells

Cellxgene allows a user to select, subset and filter cells based on gene expression cutoffs and categorical metadata attributes. There are multiple ways to do this. Here we will subset a population of 'B- cells' by different ways described below -

1. Lasso selection - Select the lasso tool and draw a lasso around the B-cell cluster. Then click on the subset icon to create the subset based on lasso selection. All other cells will disappear from the embedding plot. To bring back all the cells, click on full daatset icon next to subset icon.

![Highlight](../img/Cellxgene/15.png) 

![Highlight](../img/Cellxgene/16.png) <center>**Figure 15.** Lasso tool</center>

2. Categorical selection - Select cell type of interest from the categorical metadata on the left hand side bar and click on the subset icon.

![Highlight](../img/Cellxgene/17.png) <center>**Figure 16.** Categorical selection</center>

## Comparing the expression of multiple genes

Cellxgene explorer allows users to compare the expression of multiple genes via bivariate plots. Here, for example, we have chosed two genes - CD8A (specific to T cells) and CD14 (specific to monocytes/macrophages). First, search the genes and create a subset of associated cell types.

![Highlight](../img/Cellxgene/18.png) <center>**Figure 17.** Subset of cells selected via gene marker</center>

Then, using clipping tool, remove the outliers. When users color the embedding plot, gene expression would be evident in the associated cell type.

![Highlight](../img/Cellxgene/19.png) 

![Highlight](../img/Cellxgene/20.png) <center>**Figure 18.** Embedding plot highlighting the cell type via gene marker</center>

Users can create a bivariate plot by plotting one gene on x-axis and another on the y-axis.

![Highlight](../img/Cellxgene/21.png) <center>**Figure 19.** Bivariate plot</center>

## Finding marker genes

Cellxgene can be used to find marker genes for a cell type. Here we will compare the expression of genes in T cells versus B cells. Select T and B cells on the categorical metadata - cell type and create a subset of these populations.

![Highlight](../img/Cellxgene/22.png) <center>**Figure 20.** Subset of selected cell types</center>

Using lasso tool, select the cells and set the population 1 on the toolbar. Similarly define the second population.

Note: Cellxgene can only analyze 50000 cells so each population should not exceed this number.

![Highlight](../img/Cellxgene/23.png) <center>**Figure 21.** Differential expression icon</center>

Click on differential expression icon in the toolbar and a list of genes differentially expressed will appear in the right-hand side bar.

![Highlight](../img/Cellxgene/25.png) <center>**Figure 22.** Differential expression between two different cell types</center>

Click on the droplet icon to observe the gene expression in the two populations of cells.

![Highlight](../img/Cellxgene/26.png) 

![Highlight](../img/Cellxgene/27.png) <center>**Figure 23.** Cell populaion highlighted by differentially expressed genes</center>

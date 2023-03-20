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

![](RackMultipart20230318-1-tfv7xc_html_4c8fee45bc87872b.png)

![](RackMultipart20230318-1-tfv7xc_html_92048e632e84f32d.png)

## Examining categorical metadata

Categorical data can be examined in multiple ways -

1. Colour embedding plots - The plot can be colored by different metadata. The plot can be colored by clicking on the drop icon next to metadata of interest. Color code is displayed next to the cell count of the chosen metadata.
2. Cell counts - Each metadata category will have the cell counts along with the color code.

![](RackMultipart20230318-1-tfv7xc_html_6e4f951546772a4e.png)

Example 1 - Embedding colored by cell type

![](RackMultipart20230318-1-tfv7xc_html_f5f5578867e8ebf.png)

Example 2 - Embedding colored by disease

1. Making the selection of cells - Upon hovering over the metadata on left side, user can see the type of cells highlighted in the embedding plot.

![](RackMultipart20230318-1-tfv7xc_html_1097b287dbb99187.png) Users can click on the 'display categories' icon to display the labels over the cluster centroid.

![](RackMultipart20230318-1-tfv7xc_html_e9033ddf6e8bb649.png)

Users can select/deselect the cells by choosing the checkbox next to the value in the categorical metadat field. The unselected cells will have a smaller point size on the embedding plot.

![](RackMultipart20230318-1-tfv7xc_html_cac9c638004bf8aa.png)

Users can use the checkmark on the parent metadata to delect all the cells and then select the cells of interest, this makes it easier to highlight the cells that pertain to a specific tissue.

![](RackMultipart20230318-1-tfv7xc_html_f7c290c8a64b9f8.png)

1. Viewing the relationship between different categorical fields - After users color by a particular metadata category, for example 'Cell type', users can see the distribution of the cell type in any other category by expanding that categorical metadata field, for example 'tissue'. In the plot below, cell types belonging to a particular tissue are clearly shown with the colors in the bar.

![](RackMultipart20230318-1-tfv7xc_html_976292bf023d8448.png)

## Finding cells where gene is expressed

Numerical metadata on the embedding plot can be examined and used to filter and select cells. The numerical data is present on the right hand sidebar and users can click on the droplet icon to color the plot by qc metrics (for example, n\_genes\_by\_counts - number of genes that have been detected in the cell).

![](RackMultipart20230318-1-tfv7xc_html_177cae3eca678924.png)

To deal with outliers, users can clip the data by clicking on 'clip' icon on top of the toolbar. For example, here we have set the values to 0 to 99 percentile and we can observe the change in the scale, graph and the embedding plot.

![](RackMultipart20230318-1-tfv7xc_html_d1ca9b96faa2ade8.png)

To find the gene expression pattern for a particular gene, type the gene name in the search bar.

![](RackMultipart20230318-1-tfv7xc_html_a2443f582b31deb9.png)

For example, if we search for STAT3, a graph will appear on the top right side that contains

- Gene name
- A histogram depicting the gene expression
- Remove icon
- X-plot/Y-plot for bivariate plots
- Droplet icon to color the cells based on the gene expression

![](RackMultipart20230318-1-tfv7xc_html_462af765b244cf9a.png)

Users can use the clipping tool to remove the outliers to understand the gene expression better.

![](RackMultipart20230318-1-tfv7xc_html_44909d1edfde7c30.png)

## Selecting and subsetting cells

Cellxgene allows a user to select, subset and filter cells based on gene expression cutoffs and categorical metadata attributes. There are multiple ways to do this. Here we will subset a population of 'B- cells' by different ways described below -

1. Lasso selection - Select the lasso tool and draw a lasso around the B-cell cluster. Then click on the subset icon to create the subset based on lasso selection. All other cells will disappear from the embedding plot. To bring back all the cells, click on full daatset icon next to subset icon.

![](RackMultipart20230318-1-tfv7xc_html_d0901f208c9c8520.png)

![](RackMultipart20230318-1-tfv7xc_html_bfb9be3814abaa10.png)

1. Categorical selection - Select cell type of interest from the categorical metadata on the left hand side bar and click on the subset icon.

![](RackMultipart20230318-1-tfv7xc_html_708fbec566af7c1a.png)

## Comparing the expression of multiple genes

Cellxgene explorer allows users to compare the expression of multiple genes via bivariate plots. Here, for example, we have chosed two genes - CD8A (specific to T cells) and CD14 (specific to monocytes/macrophages). First, search the genes and create a subset of associated cell types.

![](RackMultipart20230318-1-tfv7xc_html_a6e0bd275767df6f.png)

Then, using clipping tool, remove the outliers. When users color the embedding plot, gene expression would be evident in the associated cell type.

![](RackMultipart20230318-1-tfv7xc_html_3783318e7aa04eba.png)

![](RackMultipart20230318-1-tfv7xc_html_6eca7222e6f34841.png)

Users can create a bivariate plot by plotting one gene on x-axis and another on the y-axis.

![](RackMultipart20230318-1-tfv7xc_html_75774f537f7901c9.png)

## Finding marker genes

Cellxgene can be used to find marker genes for a cell type. Here we will compare the expression of genes in T cells versus B cells. Select T and B cells on the categorical metadata - cell type and create a subset of these populations.

![](RackMultipart20230318-1-tfv7xc_html_40741c5b234b177c.png)

Using lasso tool, select the cells and set the population 1 on the toolbar. Similarly define the second population.

Note: Cellxgene can only analyze 50000 cells so each population should not exceed this number.

![](RackMultipart20230318-1-tfv7xc_html_41cb367fed645743.png)

Click on differential expression icon in the toolbar and a list of genes differentially expressed will appear in the right-hand side bar.

![](RackMultipart20230318-1-tfv7xc_html_a0ffe4c1996ea849.png)

![](RackMultipart20230318-1-tfv7xc_html_f68aa9f141ffe8b9.png)

Click on the droplet icon to observe the gene expression in the two populations of cells.

![](RackMultipart20230318-1-tfv7xc_html_4656ba7df37ab547.png)

![](RackMultipart20230318-1-tfv7xc_html_d3b3ab9ccdc36449.png)

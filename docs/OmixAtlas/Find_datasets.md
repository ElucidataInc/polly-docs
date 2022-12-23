
# 1. How to find datasets?

## 1.1 Search based on operator hints

![Search Bar](../img/OmixAtlas-Images/1a.png) <center>**Figure 1.** Search Bar</center>

This search bar is driven by Elasticsearch where users can search for keywords which are present across:

- Source metadata fields such as title, description and overall design and
- Curated metadata fields such as tissue, drug, cell line, cell type, disease, organism, gene, gene modification and dataset id.

It allows fuzzy search as well. For example, "transcriptomics" in the search keyword will show results for "transcriptome" or "transcript" as well.

This search bar supports the following operations to help users with some advanced operations such as AND (&), OR (|), NOT (~), EXACT ("text")

While parsing the search keyword, the algorithm assigns the following priority for different operators:- Brackets () \> AND (&) \> OR (|) \> NOT (~)

- Use brackets ( ) to ensure the operators in your query get executed in the exact order you want.
- term1 & term2 | term3 can be executed in two ways:

(term1 & term2) | term3 or term1 & (term2 | term3)

Using one of the above two will help remove ambiguity.

● Use exact matches " " to do stricter spelling matches with search keywords

○ Eg. transcriptome vs "transcriptome"

The former will match the word "transcriptome" & also words closer to

it such as "transcriptomics" or "transcript", whereas

The latter would only match the exactly spelled word "transcriptome"

○ Eg. renal cancer vs "renal cancer"

The former will behave same as (renal | cancer) with a slight difference:

● (renal | cancer) can match a dataset containing "renal" & "cancer" in any single field or even different fields

● renal cancer will return those datasets at top which contain both "renal" & "cancer" in the same field, followed by the ones which contain only one of these words in any single field.

The latter would only match a dataset with a field containing exact word

"renal" followed by the exact word "cancer", in the same order.

● Use operator NOT ~ to exclude unwanted results from your search

○ Eg. (HCC & Hepatocellular carcinoma) ~ "radiotherapy"

The above search would match datasets which contains both "hcc" and at least one of the words ("hepatocellular"or "carcinoma").

Out of those, ignore the datasets containing the word "radiotherapy".

Some of the Search examples are as follows:-

- (ITIH1 upregulation | FN1 downregulation) & "Fatty liver"

![Example](../img/OmixAtlas-Images/2a.png) <center>**Figure 2.** Example 1</center>

- (Hepatocellular carcinoma | HCC) ~ Radiotherapy

![Example](../img/OmixAtlas-Images/3a.png) <center>**Figure 3.** Example 2</center>


- CDK7 & "CBM signaling pathway" & inhibition

![Example](../img/OmixAtlas-Images/4a.png) <center>**Figure 4.** Example 3</center>


- (somatic mutation) & (hepatocellular | renal) & (cancer | carcinoma)

![Example](../img/OmixAtlas-Images/5a.png) <center>**Figure 5.** Example 4</center>

## 1.2 Search and Filter based on Ontology

After shortlisting the datasets using a search bar, users can further find the desired datasets using the Filter Datasets function on the left. The total number of datasets related to a filter can be seen in the brackets next to it. These filters are configurable and may be different for different OmixAtlasses.

![Filter](../img/OmixAtlas-Images/6a.png) <center>**Figure 6.** Filtering Datasets</center>

For a standard OmixAtlas, following filters are available

- Cell Line

Users can filter datasets by searching for desired cell lines using following three options from a drop down menu -

- Cell lines
- Disease cell line
- Tissue specific cell line

If users choose disease, cell lines related to those diseases are shown and not the disease itself.

![Filter](../img/OmixAtlas-Images/7a.png) <center>**Figure 7.** Filter by Cell Line</center>


- Disease

Users can filter the datasets by selecting the disease of interest.

![Filter](../img/OmixAtlas-Images/8a.png) <center>**Figure 8.** Filter by Disease</center>


- Tissue

Users can filter the datasets by selecting the tissue of interest.

![Filter](../img/OmixAtlas-Images/9a.png) <center>**Figure 9.** Filter by Tissue</center>

- Organism

Users can filter the datasets by selecting the organism of interest.

![Filter](../img/OmixAtlas-Images/10a.png) <center>**Figure 10.** Filter by Organism</center>

- Data Type

Users can filter the datasets by selecting the data type.

![Filter](../img/OmixAtlas-Images/11a.png) <center>**Figure 11.** Filter by Data Type</center>

- Drug

Users can filter the datasets by selecting the drug of interest.

![Filter](../img/OmixAtlas-Images/12a.png) <center>**Figure 12.** Filter by Drug</center>

- Platform

Users can filter the datasets by selecting the platform.

![Filter](../img/OmixAtlas-Images/13a.png) <center>**Figure 13.** Filter by Platform</center>


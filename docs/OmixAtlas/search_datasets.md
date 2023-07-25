
# 1. How to find datasets?

## 1.1 Search based on operator hints

![Search Bar](../img/OmixAtlas-Images/1a.png) 

Elasticsearch powers this search engine, allowing users to look for keywords across:

-   Source metadata fields such as **title, description,** and **overall design**

-   Curated metadata fields such as **tissue, drug, cell line, cell type, disease, organism, gene, gene modification,** and **dataset id**

Some of the features of the search are - 

-   It allows fuzzy search as well. For example, "transcriptomics" in the search keyword will show results for "transcriptome" or "transcript" as well.

-   This search bar supports the following operations to help users with some advanced operations such as AND (&), OR (|), NOT (~), EXACT ("text")

### Use of operators
-   While parsing the search keyword, the algorithm assigns the following priority to different operators:- **Brackets () > AND (&) > OR (|) > NOT (~)**
<ul>
<li>Use brackets ( ) to ensure the operators in your query get executed in the exact order you want.
<ul>
   <li> term1 & term2 | term3 can be executed in two ways: (term1 & term2) | term3 or term1 & (term2 | term3). 
    <li>Using one of the two above will help remove ambiguity.
    </li>
</ul>
<li>Use exact matches " " to do stricter spelling matches with search keywords
  <ul>
    <li>Eg. transcriptome vs "transcriptome" - The former will match the word "transcriptome" & also words closer to it, such as "transcriptomics" or "transcript", whereas the latter will only match the exactly spelled word "transcriptome".
<li>Eg. **renal cancer vs "renal cancer"** - The former will behave the same as (renal | cancer) with a slight difference:
<ul><li>(renal | cancer) can match a dataset containing "renal" & "cancer" in any single field or even different fields.
<li>renal cancer will return those datasets at the top which contain both "renal" & "cancer" in the same field, followed by those that contain only one of these words in any single field.
  <li>The latter would only match a dataset with a field containing the exact word "renal" followed by the exact word "cancer", in the same order.
</li></li></ul>
    </li>
  </li>
</ul>
</li>
</ul>
  <ul><li>Use operator NOT ~ to exclude unwanted results from your search
  <ul><li>Eg. (HCC & Hepatocellular carcinoma) ~ "radiotherapy"
<li>The above search would match datasets that contain both "hcc" and at least one of the words ("hepatocellular"or "carcinoma").
<li>Out of those, ignore the datasets containing the word "radiotherapy".</li></li></ul>

### Search Examples

Some of the Search examples are as follows:-

- (ITIH1 upregulation | FN1 downregulation) & "Fatty liver"

![Example](../img/OmixAtlas-Images/2a.png) 

- (Hepatocellular carcinoma | HCC) ~ Radiotherapy

![Example](../img/OmixAtlas-Images/3a.png) 


- CDK7 & "CBM signaling pathway" & inhibition

![Example](../img/OmixAtlas-Images/4a.png) 


- (somatic mutation) & (hepatocellular | renal) & (cancer | carcinoma)

![Example](../img/OmixAtlas-Images/5a.png) 

### Search based on Pubmed ID

All datasets include the PubMed ID for easy access to the publication. Additionally, the search bar allows users to look for datasets based on the PubMed ID.

![OA_search_5](https://github.com/ElucidataInc/polly-docs/assets/107244183/88130d99-4bb1-42c8-831e-62d0eabd9848)

## OmixAtlas Searching and Filtering 
[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/bHdl6I2YAoc/0.jpg)](https://www.youtube.com/watch?v=bHdl6I2YAoc)


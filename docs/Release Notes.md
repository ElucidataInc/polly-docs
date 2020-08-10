#Release Notes

<!--July 17th, 2020-->

<details open>
<summary><font size="+1"><b>July 17th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>We have released the newest version of <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">Polly CLI v0.1.18</a> enabling you to run a CLI job without the need of "secret" key if the private docker is on Polly.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li><a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow.</a>has N and C as indistinguishable isotopes</li>
    <li>Improved the stability of both Shiny and Desktop Appications.</li>
    <li>Communication within the infrastructure is now through encrypted keys.</li>
    <li>Shiny apps as well as shiny states are encrypted during transit as well as storage.</li>
    <li>Added encryption for the disks running the computations.</li>
    <li>Encrpted buckets containing credentials.</li>
  </ul>
</details>

<hr>

<!--July 3rd, 2020-->
<details>
<summary><font size="+1"><b>July 3rd, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
    <li>We have released the newest version of <a href="https://github.com/ElucidataInc/ElMaven/releases/tag/v0.11.0">El-MAVEN v0.11.0.</a></li>
    <li>Polly now provides its own docker repository for easy <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html#docker-management">management of dockers</a>.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Added Si as an indistinguishable isotope in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html">Labeled LC-MS Workflow.</a></li>
    <li>Introduced pre-processing functionlaities along with updated selections and heatmap for visualization in <a href="https://docs.elucidata.io/Apps/Lipidomics%20Data/Lipidomics%20Visualization%20Dashboard.html">Lipidomics Visualziation Dashboard</a>.</li>
  </ul>
  <p class="Deprecated-button">Deprecated</p>
  <ul>
    <li>Deprecated El-MAVEN FirstView Integration.</li>
  </ul>
</details>

<hr>

<!--June 19th, 2020-->
<details>
<summary><font size="+1"><b>June 19th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
   <li>We now support reactions from <i>Drosophila melanogaster</i> for integrated pathway analysis in IntOmix.</li>
    <li>Introduced <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html#statistical-analysis">pathway enrichment and pathway view</a> feature along with <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Dual%20Mode%20Visualization.html#comparative-analysis">comparative analysis</a> in Dual Mode Data Visualization.</li>
    <li>DEPMAP CCLE (DEPMAP Cancer cell line expression data and dependency scores for genes) repository has been added in <a href="https://docs.elucidata.io/Data%20Lake.html">Data Lake</a>.</li>
    <li>Implemented input file access from the sub-folders of a project for applications.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>The Single Cell Downstream docker is updated with these new packages: rpy2, anndata2ri (Python packages), ExperimentHub (R package).</li>
    <li>Added a GPU instance for <a href="https://docs.elucidata.io/Scaling%20compute/Polly%20CLI.html">Polly CLI</a>.</li>
  </ul> 
</details>

<hr>

<!--June 5th, 2020-->
<details>
<summary><font size="+1"><b>June 5th, 2020</b></font></summary>
<br>
  <p class="new-button">New</p>
  <ul>
   <li>Introduced visualization of labels in <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MS%20Workflow.html#visualization">stacked plot</a> within Labeled LC-MS Workflow.</li>
    <li>Enabled least privilege access for stringent access policies.</li>
    <li>Encryption of data in transit and at rest.</li>
  </ul>
  <p class="update-button">Update</p>
  <ul>
    <li>Improved access logs throughout the platform.</li>
    <li>Enhanced security using a secrets management service.</li>
    <li>Implemented regular backups and versioning of data.</li>
  </ul> 
</details>

<hr>

<!--May 22nd, 2020-->
<details>
  <summary><font size="+1"><b>May 22nd, 2020</b></font></summary>
  <br>
  <p class="new-button">New</p>
  <ul>
    <li>Introduced <a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/CompoundDiscoverer%20QuantFit.html">Polly QuantFit</a> node in <a href="https://www.thermofisher.com/in/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/compound-discoverer-software.html">Compound Discoverer<sup>TM</sup></a> that allows peak picking and absolute quantification on raw data obtained from a Thermo Scientific<sup>TM</sup> Mass Spec instrument.</li>
  </ul>
</details>

<hr>

<!--May 8th, 2020-->
<details>
  <summary><font size="+1"><b>May 8th, 2020</b></font></summary>
  <br>
  <p class="new-button">New</p>
  <ul>
    <li>We now host our desktop application, <a href="https://docs.elucidata.io/Apps/Metabolomic Data/El-MAVEN.html">El-MAVEN on Polly</a>.</li>
    <li><a href="https://docs.elucidata.io/Apps/Metabolomic%20Data/Labeled%20LC-MSMS%20Workflow.html#phibeta-tab">Phi calculation</a> feature has been added to Labeled LC-MS/MS Workflow.</li>
  </ul> 
  <p class="update-button">Update</p>
  <ul>
    <li>Changed the optimized color palette in IntOmix from a red-yellow-green scale to a more intuitive red-green scale. All upregulated metabolites or genes are represented by a shade of red and downregulated metabolites or genes as a shade of green.</li>
    <li>Changed the non-optimized color palette in IntOmix from a pink-purple scale to a red-green scale to remove ambiguity.</li>
  </ul> 
</details>  

<hr>

<!--April 24th, 2020-->
<details>
  <summary><font size="+1"><b>April 24th, 2020</b></font></summary>
  <br>
  <p class="new-button">New</p>
  <ul>
    <li>COVID-19 (Transcriptional datasets for SARS viruses, viral infections, and therapeutics for novel coronavirus) repository has been added in <a href="https://docs.elucidata.io/Data%20Lake.html">Data Lake</a>.</li>
    <br>
    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
    <iframe src="https://www.youtube.com/embed/AYgAb5Lbj4g" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
  </ul> 
</details>  

<hr>

<br />

<!--button style-->
<style>
  .update-button {
    background-color: #4C61AF;
    border: 1px solid #364574;
    border-radius: 70px;
    color: #FFFFFF;
    padding: 0px 5px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 12px;
    margin: 4px 2px;
    cursor: default;
  }
  .new-button {
    background-color: #4CAF50;
    border: 1px solid #367437;
    border-radius: 70px;
    color: #FFFFFF;
    padding: 0px 5px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 12px;
    margin: 4px 2px;
    cursor: default;
  }
  .Deprecated-button {
    background-color: #b30000;
    border: 1px solid #b30000;
    border-radius: 70px;
    color: #FFFFFF;
    padding: 0px 5px;
    text-align: center;
    text-decoration: none;
    display: inline-block;
    font-size: 12px;
    margin: 4px 2px;
    cursor: default;
  }
</style>

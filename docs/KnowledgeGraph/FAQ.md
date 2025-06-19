# FAQs & Troubleshooting

## Common Issues & Solutions

### 1. I can’t see the application in the Applications tab on Polly. What should I do?
Ensure that your email address is mapped to the application. If you're unsure, please contact Polly Support at [polly.support@elucidata.io](mailto:polly.support@elucidata.io) for access.

### 2. The application is taking too long to load. How can I fix this?
Try closing and relaunching the application. If the problem continues, contact Polly Support.

### 3. The graph is not generating. What could be wrong?
Make sure you have selected at least two nodes before clicking **Generate**. It's also possible that the selected nodes are not connected.

### 4. The graph looks cluttered. How can I make it easier to view?
Use the zoom controls or apply filters to adjust the view for better readability.

### 5. I’m unable to download the data. What should I check?
Ensure that your browser allows pop-ups and downloads. You may also need to disable ad blockers for Polly and its subdomains.

---

## Frequently Asked Questions (FAQs)

### 1. Which resources are available in the Knowledge Graph (KG)?
Refer to the _KG Data and Supported Resources_ section of this document for a complete list.

### 2. Can I access KG data programmatically?
Yes. KG data can be accessed in read-only mode using the `polly-python` package, either locally or within Polly Notebooks. More details are available [here](https://docs.polly.elucidata.io/polly-python/PollyKG.html).

### 3. I don’t see my entity of interest or a specific category in the dropdown. Why?
This could be due to the entity or category not being available in the KG data.

Refer to the _Selecting Nodes & Generating a Knowledge graph section_ for more details.

### 4. Is there a limit to how many hops are shown in the graph?
Yes. Currently, the graph displays a maximum of **1-hop** relationships.

### 5. How can I edit the Cypher query in Polly Co-Scientist?
Currently, Polly Co-Scientist does not support direct in-place editing of Cypher queries. However, User can copy the generated Cypher query from the output, make desired edits, paste the modified query into the input box, and submit it. Polly will recognize it as a Cypher query and execute it.

### 6. How can I download the data generated in Polly Co-Scientist?
You can download the results of your query in two formats: **CSV** and **JSON**.  
- Simply use the download option provided on the right side of the output section to get the data.

### 7. Can I continue my search from a previous query or session?  
Currently, Polly Co-Scientist does not support session memory or query chaining. Each query is treated as a standalone input, so:
- Previous queries and their results are **not retained**, and
- New queries cannot reference earlier ones.

### 8. Are there any browser or hardware requirements for using this application?
Yes. For the best experience:

- Use the latest version of **Chrome** or **Firefox**. Some features may not function correctly in other browsers.
- Disable ad blockers on Polly and its subdomains.
- Use a system with at least **8 GB RAM** and a processor equivalent to or better than **Intel Core i5**.

---

## Getting Help

For additional support, contact the Polly KG Explorer team at [polly.support@elucidata.io](mailto:polly.support@elucidata.io).
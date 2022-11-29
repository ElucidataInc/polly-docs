class Curation:
    """
The Curation class contains wrapper functions around the models used for 
semantic annotations of string/text.

Curation functions are able to recognise different entities given a text, normalise them based on certain nomenclature such as Polly compatible ontologies. Entities that  are supported are:  "disease", "drug", "species", "tissue", "cell_type", "cell_line", "gene".

Args:
   token (str): token copy from polly.
    
Usage:
	from polly.curation import Curation

 	curationObj = Curation(token)
    
    """
 
    def standardise_entity(
        self,
        mention: str,
        entity_type: str,
        context: Optional[str] = None,
        threshold: Optional[float] = None,
    ) -> dict:
        """
	Map a given mention (keyword) to an ontology term.

	Given a text and the type of entity it is, users can get the Polly compatible ontology for the text such as the MESH ontology. 
        
	Args:
            mention (str): mention of an entity e.g. "Cadiac arrythmia"
            entity_type (str): Should be one of
            ['disease', 'drug', 'tissue', 'cell_type', 'cell_line', 'species', 'gene']
            context (str): The text where the mention occurs.
            This is used to resolve abbreviations
        
	Returns:
            dict : Dictionary containing keys and values of the entity type, ontology (such as NCBI, MeSH), ontology ID (such as the MeSH ID), the score (confidence score), and synonyms if any

	Raises:
		requestException : Invalid Request
        """

    def recognise_entity(
        self,
        text: str,
        threshold: Optional[float] = None,
        normalize_output: bool = False,
    ):
        """
	Run an NER model on the given text. The returned value is a list of entities along with span info.

	Users can simply recognise entities in a given text without any ontology standardisation (unlike the annotate_with_ontology function which normalises as well).
        
	Args:
		text (str): input text
		normalize_output (bool): whether to normalize the keywords
        
	Returns:
            entities (List[dict]): returns a list of spans containing the keyword, start and end index of the keyword and the entity type

	Raises:
		requestException : Invalid Request
        """
       

    def annotate_with_ontology(
        self,
        text: str,
    ) -> List[Tuples]:

        """
        Tag a given piece of text. A "tag" is just an ontology term. Annotates with Polly supported ontologies.
        
	This function calls recognise_entity followed by normalize.
	Given a text, users can identify and tag entities in a text. Each entity/tag recognised in the text contains the 	name( word in the text identified), entity_type and the ontology_id.

        Args:
            text (str): Input text

        Returns:
            set of unique tags
        """

    def find_abbreviations(self, text: str) -> Dict[str, str]:
        """
        To run abbreviation detection separately.
        Internally calls a normaliser.
        
	Args:
            text (str): The string to detect abbreviations in
        
	Returns:
            Dictionary with abbreviation as key and full form as value

	Raises:
		requestException : Invalid Request
        """
        

    def assign_control_pert_labels(
        self, sample_metadata, columns_to_exclude=None
    ) -> pandas.DataFrame:
        """Returns the sample metadata dataframe with 2 additional columns.
	is_control - whether the sample is a control sample
	control_prob - the probability that the sample is control

        Args:
		sample_metadata (pandas.DataFrame): Metadata table
		columns_to_exclude (Set[str]): Any columns which don't play any role in determining the label, e.g. any arbitrary sample identifier

        Returns:
            DataFrame with input data frame with 2 additional columns

	Raises:
		requestException : Invalid Request
        """
    
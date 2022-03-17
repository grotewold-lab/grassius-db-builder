
class HmmscanResult:
    """
    An object containing results from running hmmscan
    
    Attributes:
    -----------
    data : DataFrame
        parsed data where each row is a matched domain-sequence pair
        includes "accession" (domain name from hmm file)
        includes "query name" (sequence header from fasta)
    footer : str
        raw text from end of blast output (starting at the first line starting with "#")
    """
    
    def __init__(self, data, footer):
        """
        Constructor for use within run_hmmscan
        """
        self.data = data
        self.footer = footer
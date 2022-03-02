
class BlastResult:
    """
    An object containing results from running TBLASTN where...
    
    - the query was one protein sequence
    - the database is a fasta file containing DNA sequences
    
    Get an instance using run_tblastn() or read_blast_output()
    
    Attributes:
    -----------
    header : str
        raw text from beginning of blast output (up to the first line starting with ">")
    data : DataFrame
        parsed data where each row is a hit
        includes "Chrom","Start_Pos","Stop_Pos","Identities%",... for each hit
    footer : 
        raw text from end of blast output (starting at the first line starting with "lamdba")
    """
    
    def __init__(self, header,data,footer):
        """
        Constructor for use within read_blast_output
        """
        self.header = header
        self.data = data
        self.footer = footer
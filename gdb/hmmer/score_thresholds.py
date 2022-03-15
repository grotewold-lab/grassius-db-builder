# this file handles the process of discounting low-scoring hmmscan results

from .hmmscan_result import HmmscanResult
import pandas as pd


all_thresholds = {

    
    # test
    "PF07716":      1.0,
    
    
    # thresholds below are from iTAK 
    # https://github.com/kentnf/iTAK/blob/master/database/GA_table.txt
    "PF00096":      7.80,
    "PF00642":      10.00,
    "PF01422":      15.90,
    "VARL":         45.0,
    "Alfin-like":   88.35,
    "VOZ":          123.3,
    "ULT":          52.45,
    "HRT":          22.1,
    "STER_AP":      110.9,
    "DNC":          10.4,
    "LUFS":         38.8,
    "G2-like":      24.35,
    "Trihelix":     26.4,
    "NF-YC":        52.7,
    "NF-YB":        61.0,
    "STAT":         150.0,
    "WUS-HB":       45.2
}


def get_filtered_hmmscan_result( hmmscan_result ):
    """
    Get a filtered version of the given hmmscan result.
    
    returns a new instance of HmmscanResult
    
    Arguments:
    ----------
    hmmscan_result -- an instance of HmmscanResult
    """
    
    in_df = hmmscan_result.data
    out_df = pd.DataFrame( columns=in_df.columns )
    
    # start copying rows from input to output
    for row in in_df.index:
        acc_name, score = in_df.loc[row,["accession","fs_score"]]
        
        # skip low-scoring rows
        if acc_name in all_thresholds.keys():
            threshold = all_thresholds[acc_name]
            if score < threshold:
                continue
                
        out_df.loc[row,:] = in_df.loc[row,:]
    
    # return filtered data
    return HmmscanResult( out_df, hmmscan_result.footer )
    



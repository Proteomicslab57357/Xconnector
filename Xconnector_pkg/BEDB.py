import urllib
import urllib.request
import pandas as pd
from pandas import DataFrame
import numpy as np 


def GetInfo(accessions, identifier = "hmdb"):

    """

    Retrieve information about small molecule metabolites found in The Blood Exposome Database

    Args:

        accessions (list): list of HMDB ID 
        identifier (string): HMDB, CID, KEGG, formula, SMILES, or InChIKey

    Returns:

        A pandas object for Data frame(s) contains general information retrieved from The Blood Exposome Database about metabolites

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        df = GetInfo(["HMDB0000001","HMDB0000002","HMDB0000005"])
        print (df)

    """

    if type(accessions) != list:
        raise TypeError ("GetInfo takes accessions as a list as an argument")

    df_db = pd.read_csv("DB//BloodExpsomeDatabase_version_1.0.csv")

    if identifier == "HMDB":
        df_db["HMDB_ID"]= df_db["HMDB_ID"].apply(str)

        df_db = df_db.loc[df_db["HMDB_ID"].isin(accessions)]
        return(df_db)

    elif identifier == "CID":
        df_db["PubChem_CID"]= df_db["PubChem_CID"].apply(str)
        df_db = df_db.loc[df_db["PubChem_CID"].isin(accessions)]

        return(df_db)

    elif identifier == "KEGG":
        df_db["KEGG_ID"]= df_db["KEGG_ID"].apply(str)

        df_db = df_db.loc[df_db["KEGG_ID"].isin(accessions)]
        return(df_db)

    elif identifier == "formula":
        df_db["Molecular_Formula"]= df_db["Molecular_Formula"].apply(str)

        df_db = df_db.loc[df_db["Molecular_Formula"].isin(accessions)]
        return(df_db)
    
    elif identifier == "SMILES":
        df_db["CanonicalSMILES"]= df_db["CanonicalSMILES"].apply(str)

        df_db = df_db.loc[df_db["CanonicalSMILES"].isin(accessions)]
        return(df_db)
        
    elif identifier == "InChIKey":
        df_db["InChIKey"]= df_db["InChIKey"].apply(str)
        df_db = df_db.loc[df_db["InChIKey"].isin(accessions)]
        return(df_db)




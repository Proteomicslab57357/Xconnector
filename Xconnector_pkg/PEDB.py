import urllib
import urllib.request
import pandas as pd
from pandas import DataFrame
import numpy as np 


def GetInfo(accessions, identifier):

    """

    Retrieve information about molecules found in The PolyphenolExplorer DB

    Args:

        accessions (list): pubchem_compound_id, name, formula
        identifier (string): pubchem_compound_id, formula,molecules name 

    Returns:

        A pandas object for Data frame(s) contains general information retrieved from The Blood Exposome Database about metabolites

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        df = GetInfo(["C28H33O16"])
        print (df)

    """

    if type(accessions) != list:
        raise TypeError ("GetInfo takes accessions as a list as an argument")

    df_db = pd.read_csv("DB//Polyphenol_Metabolites.csv", dtype=str)

    if identifier == "pubchem_compound_id":
        #df_db["pubchem_compound_id"]= df_db["pubchem_compound_id"].apply(str)
        df_db = df_db.loc[df_db["pubchem_compound_id"].isin(accessions)]
        return(df_db)

    if identifier == "formula":
        #df_db["formula"]= df_db["formula"].apply(str)

        df_db = df_db.loc[df_db["formula"].isin(accessions)]
        return(df_db)

    if identifier == "name":
        #df_db["name"]= df_db["name"].apply(str)

        df_db = df_db.loc[df_db["name"].isin(accessions)]
        return(df_db)

def pc_getinfo(GetInfo_df):
    id_list = list(set(GetInfo_df["id"].tolist())) 
    pc_df = pd.read_csv("DB//polyphenol_classification.csv", dtype=str)

    pc_df = pc_df.loc[pc_df["compound_id"].isin(id_list)]

    return(pc_df)

def pcd_getinfo(GetInfo_df):
    id_list = list(set(GetInfo_df["id"].tolist())) 
    pcd_df = pd.read_csv("DB//Polyphenols_having_composition_data.csv", dtype=str)

    pcd_df = pcd_df.loc[pcd_df["id"].isin(id_list)]

    return(pcd_df)





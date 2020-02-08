import pandas as pd
from pandas import DataFrame
import urllib
import urllib.request
import matplotlib.pyplot as plt
import numpy as np
import re
import math

def AccData(accessions):
    """

    Retrieve MSn spectra from ReSpect database

    Args:

        accessions (list): list of ReSpect accessions

    Returns:

        Two generators object for Data frame(s) contains MSn spectra data from ReSpect database and peak data 

    Raises:

        TypeError if argument (accessions) is not a list

    Example

        for data , peak_data in AccData(["PM013507","PS058407"]):
            print (data , peak_data)

    """

    def add_all_df(dfdf_new):


        dfdf = dfdf_new.drop(["Tag"], 1)
        dfdf = dfdf.transpose()
        dfdf.columns = list(dfdf_new["Tag"])
        dfdf = dfdf.rename(index={"Data": dfdf["ACCESSION"][1]})
        dfdf = dfdf.drop(["ACCESSION"], 1 )
        new_header = []

        for i_header, y_header in zip(dfdf.iloc[0], dfdf.columns.values):
            if type(i_header) == str:
                new_header.append(f"{y_header}_{i_header}")
            else:
                new_header.append(y_header)
        
        dup = set([i for i in new_header if new_header.count(i) > 1])

        count = 0
        for n , c in enumerate(new_header):
            if c in dup:
                count += 1
                new_header[n] = f"{c}_{count}"

        dfdf.columns = new_header
        dfdf = dfdf.drop(["Sub Tag"], 0)
        return dfdf


    if type(accessions) != list:
        raise TypeError ("AccData takes list as an argument")

    main_url = "http://spectra.psc.riken.jp/menta.cgi/respect/datail/datail?"

    for i_acc in accessions:
        # creating the query search url by accessions
        dict_url = {"accession":i_acc}
        api_request =  urllib.parse.urlencode(dict_url)
        api_request = main_url + api_request
            
        # get the data in form of tables
        tables = pd.read_html(api_request)
        if len (tables) > 1 :
            tables = tables[1]
            #### modify peak data
            # save the infromation of that cell
            peak = tables.iat[-1,2]
            # remove the peak data cell values from tables
            tables.iat[-1, 2] = ""

            # creat dataframe for peak data
            peak_table = peak[peak.find("rel.int.")+9:]
            peak_table = peak_table.strip()
            peak_table = peak_table.split(" ")
            peak_table = list(filter(None, peak_table))
            peak_table = [peak_table[x:x+3] for x in range(0,len(peak_table),3)]
            peak_table[0] = [i.replace(u'\xa0', ' ') for i in peak_table[0] ]
            peak_tables = DataFrame(peak_table,columns = peak_table[0]).drop(0,0)
            peak_tables.reset_index(inplace=True, drop=True)

            rename_list_index = [i_acc] * len(peak_tables)
            index_name = list ( range(0, len(peak_tables)) )
            dict_rename_index = dict( zip( index_name , rename_list_index ) )

            peak_tables = peak_tables.rename(index=dict_rename_index)

            #remove the last row from tables
            tables.drop(index=len(tables)-1,inplace = True)
            

            tables = add_all_df(tables)
        else:
            tables = pd.DataFrame(index=[i_acc])
            peak_tables = pd.DataFrame()

        yield tables, peak_tables
      
def draw_peak(accession):
    """

    Retrieve MSn spectra from ReSpect database

    Args:

    accession (str): string of ReSpect accession

    Returns:

        BarContainer

    Raises:

        TypeError if argument (accession) is more than one or not a string

    Example

        draw_peak("PS058407").show()

    """

    if type(accession) == list and len(accession) > 1:
        raise  TypeError ("draw_peak takes only one accession")

    for _ , peak_table in AccData([accession]):

        peak_table["m/z"] = pd.to_numeric(peak_table["m/z"])
        peak_table["Relative intensity"] = pd.to_numeric(peak_table["Relative intensity"])
        peak_table.sort_values(by = "m/z", axis= 0 ,inplace=True,ascending = True)
        
    peak_table.index.name = "index"
    peak_table.reset_index(inplace=True)
    peak_table.drop("index",axis=1,inplace=True)
    peak_table.drop(["Intensity"] , axis = 1 , inplace = True)
    index_replace = dict(zip(peak_table.index , peak_table["Relative intensity"].astype(str)))
    peak_table.rename(index= index_replace, inplace = True)

    plt.style.use("seaborn-pastel")
    plt.rc('figure', figsize=(15, 15))

    peak_table.plot(x = "m/z" , y = "Relative intensity" , width=0.001,edgecolor="blue",color = "blue" ,alpha=0.6, kind = "bar" , legend=False)

    plt.suptitle('MSn spectra', fontsize=20)
    plt.xlabel('m/z', fontsize=18)
    plt.ylabel('Relative intensity', fontsize=16)
    if max(peak_table["Relative intensity"]) == 100:
        plt.ylim(0, max(peak_table["Relative intensity"])+10)
    elif max(peak_table["Relative intensity"]) == 1000:
        plt.ylim(0, max(peak_table["Relative intensity"])+100)
    else:
        plt.ylim(0, max(peak_table["Relative intensity"])+20)
    plt.title(f"ID: {accession}")

    for x, y in enumerate(peak_table["Relative intensity"]):
        plt.text(x-0.01, y + 1, str(y), color='blue', rotation=40)

    return plt

def Keyword(name,formula, exactmass,tolerance):
    """

    Searching and Retrieve MSn spectra from ReSpect database

    Args:

        name (str): compound name
        formula (str): compound formula
        exactmass (str): compound exact mass
        tolerance (str): compound tolerance

    Returns:

        Two generators object for Data frame(s) contains MSn spectra data from ReSpect database and peak data 

    Raises:

        TypeError if argument (name) is not a str
        TypeError if argument (formula) is not a str
        TypeError if argument (exactmass) is not a str
        TypeError if argument (tolerance) is not a str

    Example

        print ( Keyword("","C21H2","","") )


    """
    if type(name) != str:
        raise  TypeError ("Keyword takes str as an argument")
    if type(formula) != str:
        raise  TypeError ("Keyword takes str as an argument")
    if type(exactmass) != str:
        raise  TypeError ("Keyword takes str as an argument")
    if type(tolerance) != str:
        raise  TypeError ("Keyword takes str as an argument")


    main_url = "http://spectra.psc.riken.jp/menta.cgi/respect/search/keyword?"
    dict_url = {"name":name,"formula":formula,"exactmass":exactmass,"tolerance":tolerance,"hideGraph":"hideGraph"}
    api_request =  urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request
    xml = urllib.request.urlopen(api_request)
    contents = str(xml.read())

    acc = re.findall('menta.cgi/respect/datail/datail\?accession=[A-Z]+[0-9]+', contents)
    acc = [i[i.find("=")+1:] for i in acc]

    if len(acc) == 0:
        not_found = "Your search did not match any documents in ReSpect database"
        return not_found
    else:
        return acc



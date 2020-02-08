import urllib
import urllib.request
import pandas as pd
from pandas import DataFrame

def Geninfo(accessions):

    """

    Retrieve information about small molecule metabolites found in the human body from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains general information retrieved from LMDB about metabolites

    Raises:

        TypeError if argument (accessions) is not a list

    Note:

        If an accession has no information in LMDB the function will return an empty data frame

    Example:

        for i_data in Geninfo(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """

    def add_all_df(dfdf_new):

    
        dfdf = dfdf_new.drop(["col1"], 1)
        dfdf = dfdf.transpose()
        dfdf.columns = list(dfdf_new["col1"])
        dfdf = dfdf.rename(index={"col2": dfdf["Lmdb"][0]})

        return dfdf

    if type(accessions) != list:
        raise TypeError ("Geninfo takes list as an argument")

    main_url = "http://www.lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        # get the data in form of tables
        try:
            tables = pd.read_html(api_request)
            #removing row with no data
            tables0 = tables[0].iloc[:,0:2]
            #rename the col 
            tables0.columns = ["col1", "col2"]
  
            #take only the general information
            take_from_tables0 = ["Version","Status","Creation Date","Update Date","Lmdb","Common Name","Chemical Formula","Average Molecular Weight",
            "Monoisotopic Molecular Weight","IUPAC Name","Traditional Name","CAS Registry Number",
            "SMILES","InChI Identifier","InChI Key","Kingdom","Super Class","Class","Sub Class","Direct Parent",
            "Molecular Framework","Role","State","Cellular Locations","Biospecimen Locations","Tissue Locations",
            "DrugBank ID","FoodDB ID","Chemspider ID","KEGG Compound ID","ChEBI ID","PubChem Compound"]
            
            tables0 = tables0[tables0["col1"].isin(take_from_tables0)]
            tables0 = add_all_df(tables0)

            yield tables0

        except:
            take_from_tables0 = ["Version","Status","Creation Date","Update Date","Lmdb","Common Name","Chemical Formula","Average Molecular Weight",
            "Monoisotopic Molecular Weight","IUPAC Name","Traditional Name","CAS Registry Number",
            "SMILES","InChI Identifier","InChI Key","Kingdom","Super Class","Class","Sub Class","Direct Parent",
            "Molecular Framework","Role","State","Cellular Locations","Biospecimen Locations","Tissue Locations",
            "DrugBank ID","FoodDB ID","Chemspider ID","KEGG Compound ID","ChEBI ID","PubChem Compound"]

            not_found = ["NaN"] * len(take_from_tables0)
            not_found = dict( zip(take_from_tables0,not_found) )
            tables0 = pd.DataFrame(data=not_found , index= [0])

            tables0 = tables0.rename(index={0:i_acc })

            yield tables0


def AccData(accessions):

    """

    Retrieve information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains general information retrieved from LMDB about metabolites

    Raises:

        TypeError if argument (accessions) is not a list

    Note:

        If an accession has no information in LMDB the function will return an empty data frame

    Example:

        for i_data in AccData(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """
    def add_all_df(dfdf_new):

        dfdf = dfdf_new.drop(["col1"], 1)
        dfdf = dfdf.transpose()
        dfdf.columns = list(dfdf_new["col1"])
        dfdf = dfdf.rename(index={"col2": dfdf["Lmdb"][0]})

        return dfdf

    if type(accessions) != list:
        raise TypeError ("AccData takes list as an argument")

    main_url = "http://lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        # get the data in form of tables
        try:
            tables = pd.read_html(api_request)
        
            try:
                #removing row with no data
                tables0 = tables[0].iloc[:,0:2]
                #rename the col 
                tables0.columns = ["col1", "col2"]
                #remove thus dataframes from the big dataframe
                remove_from_tables0 = ["Synonyms","Experimental Properties","Predicted Properties","Spectra"]
                for i_remove in remove_from_tables0:
                    tables0 = tables0[tables0["col1"] != i_remove]
                #remove the Concentrations dataframe from the big dataframe
                try:
                    conc = tables0.index[tables0['col1'] == "Concentrations"].tolist()
                    tables0 = tables0[tables0["col1"] != "Concentrations"]
                    tables0.drop(index=conc[0],inplace=True)
                    tables0.drop(index=conc[0]+1,inplace=True)
                except:
                    pass
            except:
                pass
        
            tables0 = add_all_df(tables0)
            tables0 = tables0.loc[:,~tables0.columns.duplicated()]

        except:
            tables0 = pd.DataFrame(data = {"notfound":["NaN"]})
            tables0 = tables0.rename(index={0: i_acc })
            tables0.drop(["notfound"], 1 , inplace= True)

        yield  tables0


def SynonymsData(accessions):
    """

    Retrieve synonyms names information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID

    Returns:

        A generator object for Data frame(s) contains the synonyms names of metabolites from LMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in SynonymsData(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("SynonymsData takes list as an argument")

    main_url = "http://lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[1]
        except:
            tables = pd.DataFrame(data = {"Value":["NA"], "Source":["NA"]})

        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )


        yield tables


def ExpProp(accessions):
    """

    Retrieve Experimental Properties information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Experimental Properties of metabolites from LMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in ExpProp(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("ExpProp takes list as an argument")

    main_url = "http://lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[2]
        except:
            tables = pd.DataFrame(data = {"Property":["NA"], "Value":["NA"], "Reference":["NA"]    })

                
        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )
            
        yield tables

def PredProp(accessions):
    """

    Retrieve Predicted Properties information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Predicted Properties of metabolites from LMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in PredProp(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("PredProp takes list as an argument")

    main_url = "http://lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[3]
        except:
            tables = pd.DataFrame(data = {"Property":["NA"],"Value":["NA"],"Source":["NA"]})

        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )

        yield tables

def Spectra(accessions):
    """

    Retrieve Spectra information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Spectra of metabolites from LMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in Spectra(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("Spectra takes list as an argument")

    main_url = "http://lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[4].iloc[0:,0:3]
        except:
            tables = pd.DataFrame(data = {"Spectrum Type":["NA"],  "Description":["NA"], "Splash Key":["NA"] })

        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )

        yield tables


def ConcsData(accessions):
    """

    Retrieve Concentrations information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Concentrations of metabolites from LMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in ConcsData(["LMDB00001","LMDB00002","LMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("ConcsData takes list as an argument")

    main_url = "http://lmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[5].iloc[0:,0:6]
            tables.rename(columns={"Reference": "Pubmed"},inplace=True)
        except:
            tables = pd.DataFrame(data = {"Biofluid":["NA"], "Status":["NA"], "Value":["NA"], "Condition":["NA"], "Species":["NA"], "Pubmed":["NA"]})

        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )

        yield tables


def BioBro(status,biofluid,metabolite_list=list(),disease_list=list()):
    """

    Browsing Biofluids information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        status (list): list of keywords to filter with metabolite status (quantified, detected, or/and expected)

        biofluid (list): list of keywords to filter with biofluid type (other_fluids, urine, milk, plasma, serum, feces or/and ruminal_fluid)

        metabolite_list (list): default = empty list. list of LMDB ID and/or metabolite names

        disease_list (list): default = empty list. list of disases names
        
    Returns:

        A generator object for Data frame(s) contains the browsing biofluids information of metabolites from LMDB

    Raises:

        TypeError if argument (status) is not a list and out of this list [quantified, detected, or/and expected] 

        TypeError if argument (biofluid) is not a list and out of this list [other_fluids, urine, milk, plasma, serum, feces or/and ruminal_fluid] 

    Example:

        for i_result in  ( Biofluids(["quantified"],["other_fluids"],["LMDB00001","LMDB00003"]) ):
            print (i_result)

    """

    #make sure that status is fine to go
    if type (status) != list:
        raise TypeError ("status should be a list")
    opp_status= ["quantified","detected","expected"]
    if len([ i for i in status if i in opp_status ]) != len(status):
        raise TypeError ("status should be a list from ( quantified, detected or/and expected ) ") 
    
    
    #make sure that biofluid is fine to go
    if type (biofluid) != list:
        raise TypeError ("status should be a list")
    opp_biofluid= ["other_fluids","urine","milk","plasma","serum","feces","ruminal_fluid"]
    if len([ i for i in biofluid if i in opp_biofluid ]) != len(biofluid):
        raise TypeError ("biofluid should be a list from ( other_fluids,urine,milk,plasma,serum,feces and ruminal_fluid )") 
    
    
    main_url = "http://lmdb.ca/biofluids?"
    dict_url = {"utf8":"✓","filter":"true","DataTables_Table_0_length":100}

    for i_status in status:
        dict_url[i_status] = "1"
    for i_biofluid in biofluid:
        dict_url[i_biofluid] = "1"
    
    if len(metabolite_list) != 0:
        dict_url["metabolite"] = "1"
        metabolites = "\n".join(metabolite_list)
        dict_url["metabolite_list"] = metabolites

    if len(disease_list) != 0:
        dict_url["disease"] = "1"
        diseases = "\n".join(disease_list)
        dict_url["disease_list"] = diseases

    api_request =  urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request
    tables = pd.read_html(api_request)

    for i_tables in range(len(tables)):
        yield tables[i_tables] 


def structures(accession):
    """

    Retrieve Structure image of small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)


    Args:

        accessions (str): str of only one LMDB ID 

    Returns:

        png image for the metabolite structure

    Raises:

        TypeError if argument (accession) is not a string of only one LMDB ID

    Example:

        structures("LMDB00001")
   
    """
    if type(accession) != str:
        raise TypeError ("structures takes string as an argument with only one LMDB ID")

    main_url = "http://lmdb.ca/structures/"
    # creating the query search url by accessions to get the image 
    api_request = f"{main_url}{accession}/image.png" 
    image_name = f"{accession}.png"
    image= urllib.request.urlretrieve(api_request, image_name)

    return image

def MetaboBro(biofluid=list(),species=list()):
    """

    Browsing metabolites information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        biofluid (list): list of keywords to filter with biofluid type (quantified, detected, or/and expected)

        species (list): list of keywords to filter with species  (bovine, ovine, caprine, equine and/or porcine)
        
    Returns:

        Data frame contains the browsing metabolites information of metabolites from LMDB

    Raises:

        TypeError if argument (biofluid) is not a list and out of this list [other_fluids, urine, milk, plasma, serum, feces or/and ruminal_fluid] 

        TypeError if argument (species) is not a list and out of this list [bovine, ovine, caprine, equine and/or porcine] 

    Example:

        MetaboBro(biofluid=["urine"], species=["bovine"]) 

    """

    #make sure that biofluid is fine to go
    if type (biofluid) != list:
        raise TypeError ("status should be a list")
    opp_biofluid= ["other_fluids","urine","milk","plasma","serum","feces","ruminal_fluid"]
    if len([ i for i in biofluid if i in opp_biofluid ]) != len(biofluid):
        raise TypeError ("biofluid should be a list from ( other_fluids,urine,milk,plasma,serum,feces and/or ruminal_fluid )")
        
    #make sure that biofluid is fine to go
    if type (species) != list:
        raise TypeError ("status should be a list")
    opp_species= ["bovine","ovine","caprine","equine","porcine"]
    if len([ i for i in species if i in opp_species ]) != len(species):
        raise TypeError ("species should be a list from (bovine, ovine, caprine, equine and/or porcine)") 
    
    #the main url
    main_url = "http://lmdb.ca/metabolites?"
    dict_url = {"utf8":"✓","filter":"true"}

    #filters
    if len(biofluid) != 0:
        for i_bio in biofluid:
            dict_url[i_bio] = "1"

    if len(species) != 0:
        for i_species in species:
            dict_url[i_species] = "1"
    
    
    all_tables = pd.DataFrame()
    
    for page in range(1,10000):

        dict_url["page"] = page
        api_request =  urllib.parse.urlencode(dict_url)
        api_request = main_url + api_request
        tables = pd.read_html(api_request)
        tables = tables[0]
        if len(tables) != 0:
            all_tables = all_tables.append(tables, ignore_index=True, sort = False) 
        else:
            break

    return all_tables


def ClassBro(status=list(),biofluid=list()):
    """

    Browsing Classification information about small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        status (list): list of keywords to filter with metabolite status (quantified, detected, or/and expected)

        biofluid (list): list of keywords to filter with biofluid type (other_fluids, urine, milk, plasma, serum, feces or/and ruminal_fluid)
        
    Returns:

        Data frame contains the browsing Classification information of metabolites from LMDB

    Raises:

        TypeError if argument (status) is not a list and out of this list [quantified, detected, or/and expected] 

        TypeError if argument (biofluid) is not a list and out of this list [other_fluids, urine, milk, plasma, serum, feces or/and ruminal_fluid] 

    Example:

        print (ClassBro().head())
    
    """

       #make sure that status is fine to go
    if type (status) != list:
        raise TypeError ("status should be a list")
    opp_status= ["quantified","detected","expected"]
    if len([ i for i in status if i in opp_status ]) != len(status):
        raise TypeError ("status should be a list from ( quantified, detected or/and expected ) ") 
    
    
    #make sure that biofluid is fine to go
    if type (biofluid) != list:
        raise TypeError ("status should be a list")
    opp_biofluid= ["other_fluids","urine","milk","plasma","serum","feces","ruminal_fluid"]
    if len([ i for i in biofluid if i in opp_biofluid ]) != len(biofluid):
        raise TypeError ("biofluid should be a list from ( other_fluids,urine,milk,plasma,serum,feces and ruminal_fluid )") 
    

    #the main url
    main_url = "http://lmdb.ca/classyfication?"
    dict_url = {"utf8":"✓","filter":"true"}

    #filters
    if len(biofluid) != 0:
        for i_bio in biofluid:
            dict_url[i_bio] = "1"

    if len(status) != 0:
        for i_species in status:
            dict_url[i_species] = "1"

    api_request =  urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request
    tables = pd.read_html(api_request)
    tables = tables[0]

    if len(tables) == 0:
        no_data = "No Data available for this query"
        return no_data
    else:
        return tables

def ChemQuery(start=100,end=200,search_type="molecular",filters=list()):
    """

    Search by molecular weight for small molecule metabolites found in different livestock species from The Livestock Metabolome Database (LMDB)

    Args:

        start (int): default = 100, the start of molecular weight to search from.

        end (int): default = 200, the end of molecular weight to search to.

        search_type: default = "molecular", either by the Molecular weight / Average mass (molecular) or by the Monoisotopic mass (monoisotopic)

        filters (list): list of keywords to filter with metabolite status (quantified, detected, or/and expected)

    Returns:

        Data frame contains the ChemQuery search result of metabolites from LMDB

    Raises:

        TypeError if argument (start) is not a int

        TypeError if argument (end) is not a int

        TypeError if argument (search_type) is not a string 

        TypeError if argument (filters) is not a list and out of this list [quantified, detected, or/and expected] 

    Example:

        data_frame = ChemQuery(start=100,end=200,search_type="molecular",filters=["quantified"])
        print (data_frame.head())
    
    """

    #make sure that start and end is fine to go
    if type(start) != float or type(end) != float or start >= end:
        raise TypeError ("start and end arguments must be integer, and also the start must be less than the end")
    #make sure that search_type is molecular or monoisotopic
    if search_type != "molecular" and search_type != "monoisotopic":
        raise TabError ( "search_type argument must be either (molecular) for Molecular weight / Average mass, or (monoisotopic) for Monoisotopic mass" )
    #make sure that filters is fine to go
    if type (filters) != list:
        raise TypeError ("filters should be a list")
    opp_filters= ["quantified","detected","expected"]
    if len([ i for i in filters if i in opp_filters ]) != len(filters):
        raise TypeError ("filters should be a list from ( quantified, detected or/and expected ) ") 

    #the main url
    main_url = "http://lmdb.ca/structures/search/metabolites/mass?"
    dict_url = {"utf8":"✓","query_from":start,"query_to":end,"search_type":search_type,"commit":"Search"}

    #make the api url

    api_request =  urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request

    #add the filters to the api url
    if len(filters) != 0:
        for i_filter in filters:
            api_filters = f"filters[status][{i_filter}]=1"
            api_request = f"{api_request}&{api_filters}"

    def find_HMDB(xml):
        import re
        xml = str(xml)
        all_hmdb = (i.end() for i in re.finditer("LMDB",xml))
        all_hmdb = list( set ( [ f"LMDB{xml[i:i+5]}" for i in all_hmdb if xml[i:i+5].isdigit()  ] ) )
        return all_hmdb

    page = 0
    all_hmdb_id = []
    while True:
        page = page + 1
        api_request = f"{api_request}&page={page}"
        xml = urllib.request.urlopen(api_request).read()
        final_HMDB = find_HMDB(xml)
        if len(final_HMDB) != 0:
            all_hmdb_id = all_hmdb_id + [ i for i in  find_HMDB(xml)] 
        if len(final_HMDB) == 0:
            break

    return Geninfo(all_hmdb_id)


def txtsearch(query):
    """

    Advanced searching using a powerful search engine based on the Lucene query language. for more information see: http://lmdb.ca/textquery

    Args:

        query (str): a query string

    Returns:

        A generator object for Data frame(s) contains information retrieved from LMDB about metabolites

    Raises:

        TypeError if argument (query) is not a string

    Example:
    
        data = txtsearch("(histidine OR poultry) AND NOT glycolylneuraminic")
        print (data)

    """
    if type(query) != str:
        raise TypeError ("query argument should be a list")

    main_url = "http://lmdb.ca/unearth/q?"
    dict_url = {"utf8":"✓","query":query,"searcher":"metabolites"}

    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request
    xml = urllib.request.urlopen(api_request).read()
  
    def find_LMDB(xml_str):
        import re
        xml_str = str(xml_str)
        if "returned no results" not in xml_str:
            end_of_lmdb = (lmdb.end() for lmdb in re.finditer("LMDB",xml_str))
            final_LMDB = list( set( [ f"LMDB{xml_str[end:end+5]}" for end in end_of_lmdb if xml_str[end:end+5].isdigit()] ) )
            return final_LMDB
        elif "returned no results" in xml_str:
            return False

    xml_parse =  find_LMDB(xml)

    return xml_parse


def LCMS(masses,mode,adducts,tolerance,tolerance_unit):
    """

    Advanced Mass Spectrum search

    Args:

        masses (list): a query Masses (Da) maximum 700 query masses per request

        mode (str): Ion Mode (positive, negative, or neutral)

        adducts (list): Adduct Type (e.g. M+H, M+H-H2O)

        tolerance (int): Molecular Weight Tolerance ± 

        tolerance_unit (str): tolerance unit (Da or ppm)

    Returns:

        Data frame object contains information retrieved from LMDB about metabolites returned from the Mass Spectrum search

    Raises:

        TypeError if argument (masses) is not a list of length more than zero

        TypeError if argument (mode) is not a string from ("positive","negative","neutral")

        TypeError if argument (adducts) is not a list

        TypeError if argument (tolerance) is not a number

        TypeError if argument (tolerance_unit) is not a string from ("Da","ppm")

    Example:
    
        result = LCMS(masses=[300,400],mode="positive",adducts=["M+H","M+H-H2O"],tolerance=0.5,tolerance_unit="Da")
        print (result)
   
    """
    
    #make sure that masses is fine to go
    if type(masses) != list or len(masses) == 0:
        raise TypeError ("masses argument must be a list of length more than zero")
    if len(masses) > 700:
        raise TypeError ("maximum 700 query masses per request")
    #make sure that mode is fine to go
    opp_mode = ["positive","negative","neutral"]
    if type(mode) != str or mode not in opp_mode:
        raise TypeError ("mode argument should be string (positive, negative, or neutral)")
    #make sure that adduct is fine to go
    if type(adducts) != list and len(adducts) == 0:
        raise TypeError ("adducts argument must be a list")
    #make sure that tolerance is fine to go
    tolerance = float(tolerance)
    if type(tolerance) != float:
        raise TypeError ("tolerance argument must be a number")
    #make sure that tolerance_unit is fine to go
    opp_unit = ["Da","ppm"]
    if type(tolerance_unit) != str or tolerance_unit not in opp_unit:
        raise TypeError ("tolerance_unit argument should be string (Da or ppm)")

    #generate the api url
    main_url="http://lmdb.ca/spectra/ms/search?"
    dict_url= {"utf8":"✓","commit":"Search","database":"LMDB","ms_search_ion_mode":mode,"tolerance":tolerance,"tolerance_units":tolerance_unit}

    #for masses
    masses_join = "\n".join(list(map(str,masses)))
    dict_url["query_masses"] = masses_join

    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request
    #for adducts
    all_tables = []
    for api_adducts in adducts:
        api_adducts = api_adducts.replace("+","%2B")
        api_adducts = api_adducts.replace("-","%2D")
        api_request = f"{api_request}&adduct_type={api_adducts}"
        tables = pd.read_html(api_request)
        all_tables = all_tables + tables
    
    all_tables_all= pd.DataFrame()
    for i_range in range ( len(all_tables) ):

        all_tables_all = all_tables_all.append(all_tables[i_range], ignore_index=True, sort = False) 

    all_tables_all.sort_values(by=['Monoisotopic Mass'],ascending=True,inplace=True)

    return all_tables_all
    

def LCMSMS(p_ion_mass,p_ion_tolerance,parent_ion_mass_tolerance_units,ion_mode,cid,peaks,mz_tolerance,mz_tolerance_units,predicted="False"):
    """

    Advanced Tandom Mass Spectrum search

    Args:

        p_ion_mass (int): parent Ion Mass (Da)

        p_ion_tolerance (int): parent Ion Mass Tolerance ±

        ion_mode (str) : ionization Mode (positive, negative)

        parent_ion_mass_tolerance_units (str): parent Ion Mass Tolerance unit (Da or ppm)

        cid (str): CID Energy (low, med, or high)

        peaks (list): MS/MS Peak List ["M/Z Intensity"] (e.g. ["40.948 0.174","56.022 0.424"])

        mz_tolerance (int) : Mass/Charge (m/z) Tolerance ±

        mz_tolerance_units (str): Mass/Charge (m/z) Tolerance unit (Da or ppm)

        predicted  : if True it will include predicted spectra. (default="False")

    Returns:

        Data frame object contains information retrieved from LMDB about metabolites returned from Tandom Mass Spectrum search

    Raises:

        TypeError if argument (masses) is not a list of length more than zero

    Example:

        result =LCMSMS(p_ion_mass=146.0,p_ion_tolerance=0.1,parent_ion_mass_tolerance_units="Da",ion_mode="positive",cid="low",peaks=["40.948 0.174","56.022 0.424"],mz_tolerance=0.5,mz_tolerance_units="Da",predicted=True)
        print (result)

   """

    main_url = "http://lmdb.ca/spectra/ms_ms/search?"
    dict_url= {"utf8":"✓","commit":"Search","searcher":"metabolites",
    "parent_ion_mass":p_ion_mass,
    "parent_ion_mass_tolerance":p_ion_tolerance,
    "parent_ion_mass_tolerance_units":parent_ion_mass_tolerance_units,
    "ms_ms_search_ion_mode":ion_mode,
    "collision_energy_level":cid,
    "mass_charge_tolerance":mz_tolerance,
    "mass_charge_tolerance_units":mz_tolerance_units}

    if predicted == "True":
        dict_url["predicted"] = 1
    else:
        pass

    if len(peaks) != 0:
        peaks = [i.replace(" ","%20") for i in peaks]
        peaks = "%0D%0A".join(peaks)

    api_request = urllib.parse.urlencode(dict_url)
    api_request = f"{api_request}&peaks={peaks}"
    api_request = main_url + api_request

    tables = pd.read_html(api_request)
    tables = tables[0].drop(["Spectral Display Tools","Structure"],axis=1)

    return tables


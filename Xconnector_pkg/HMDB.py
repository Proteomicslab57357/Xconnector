import urllib
import urllib.request
import pandas as pd
from pandas import DataFrame



def Geninfo(accessions):

    """

    Retrieve information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID 

    Returns:

        A generator object for Data frame(s) contains general information retrieved from HMDB about metabolites

    Raises:

        TypeError if argument (accessions) is not a list

    Note:

        If an accession has no information in HMDB the function will return an empty data frame

    Example:

        for i_data in Geninfo(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """

    def add_all_df(dfdf_new):

    
        dfdf = dfdf_new.drop(["col1"], 1)
        dfdf = dfdf.transpose()
        dfdf.columns = list(dfdf_new["col1"])
        dfdf = dfdf.rename(index={"col2": dfdf["HMDB ID"][0]})

        return dfdf

    if type(accessions) != list:
        raise TypeError ("Geninfo takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
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
            take_from_tables0 = ["Version","Status","Creation Date","Update Date","HMDB ID","Common Name","Chemical Formula","Average Molecular Weight",
            "Monoisotopic Molecular Weight","IUPAC Name","Traditional Name","CAS Registry Number",
            "SMILES","InChI Identifier","InChI Key","Kingdom","Super Class","Class","Sub Class","Direct Parent",
            "Molecular Framework","Role","State","Cellular Locations","Biospecimen Locations","Tissue Locations",
            "DrugBank ID","FoodDB ID","Chemspider ID","KEGG Compound ID","ChEBI ID","PubChem Compound"]
            
            tables0 = tables0[tables0["col1"].isin(take_from_tables0)]
            tables0 = add_all_df(tables0)

            yield tables0

        except:
            take_from_tables0 = ["Version","Status","Creation Date","Update Date","HMDB ID","Common Name","Chemical Formula","Average Molecular Weight",
            "Monoisotopic Molecular Weight","IUPAC Name","Traditional Name","CAS Registry Number",
            "SMILES","InChI Identifier","InChI Key","Kingdom","Super Class","Class","Sub Class","Direct Parent",
            "Molecular Framework","Role","State","Cellular Locations","Biospecimen Locations","Tissue Locations",
            "DrugBank ID","FoodDB ID","Chemspider ID","KEGG Compound ID","ChEBI ID","PubChem Compound"]

            not_found = ["NaN"] * len(take_from_tables0)
            not_found = dict( zip(take_from_tables0,not_found) )
            tables0 = pd.DataFrame(data=not_found , index= [0])

            tables0 = tables0.rename(index={0:i_acc })

            yield tables0


def SynonymsData(accessions):

    """

    Retrieve synonyms names information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID

    Returns:

        A generator object for Data frame(s) contains the synonyms names of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in SynonymsData(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """

    if type(accessions) != list:
        raise TypeError ("SynonymsData takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
    
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

    Retrieve Experimental Properties information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID

    Returns:

        A generator object for Data frame(s) contains the  Experimental Propertiess of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in ExpProp(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("SynonymsData takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
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

    Retrieve Predicted Properties information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID

    Returns:

        A generator object for Data frame(s) contains the Predicted Properties of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in PredProp(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("SynonymsData takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
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

    Retrieve Spectra Properties information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of LMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Spectra of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in Spectra(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("Spectra takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
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


def NConcsData(accessions):
    """

    Retrieve Normal Concentrations Properties information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Normal Concentrations of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in NConcsData(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("ConcsData takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[6].iloc[0:,0:6]
            tables.rename(columns={"Reference": "Pubmed"},inplace=True)
        except:
            tables = pd.DataFrame(data = {"Biospecimen":["NA"], "Status":["NA"], "Value":["NA"], "Age":["NA"], "Sex":["NA"], "Condition":["NA"]})

        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )

        yield tables

def AConcsData(accessions):
    """

    Retrieve Abnormal Concentrations Properties information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID 

    Returns:

        A generator object for Data frame(s) contains the Abnormal Concentrations of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in AConcsData(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("ConcsData takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[7].iloc[0:,0:6]
            tables.rename(columns={"Reference": "Pubmed"},inplace=True)
        except:
            tables = pd.DataFrame(data = {"Biospecimen":["NA"], "Status":["NA"], "Value":["NA"], "Age":["NA"], "Sex":["NA"], "Condition":["NA"]})
        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )

        yield tables


def Pathways(accessions):
    """

    Retrieve Pathways names about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        accessions (list): list of HMDB ID

    Returns:

        A generator object for Data frame(s) contains the Pathways names of metabolites from HMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in Pathways(["HMDB0000001","HMDB0000002","HMDB0000005"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("SynonymsData takes list as an argument")

    main_url = "http://www.hmdb.ca/metabolites/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[5]
            for index in range(len(tables.columns)):
                try:
                    tables.drop(tables.columns[1], axis=1, inplace=True)
                except:
                    break
            tables.rename(columns={tables.columns[0]:i_acc},inplace=True)
        except:
            tables = pd.DataFrame(data = {i_acc:["NaN"]})
        yield tables


def BroDis(status=list(),biospecimen=list(),metabolite=list(),disease=list(),inborn_errors=False):
    """

    Browsing diseases information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        status (list): list of keywords to filter by metabolite status (quantified, detected, expected, or/and predicted) (default all)

        biospecimen (list): list of keywords to filter by biospecimen (blood,saliva,urine,csf,feces,sweat,breast_milk,bile,amniotic_fluid, or/and other_fluids)

        metabolite_list (list): default = empty list. list of LMDB ID and/or metabolite names

        disease_list (list): default = empty list. list of disases names

        inborn_errors (bool): default = False. if True Filter by disease type (Inborn Errors of Metabolism) is active 

    Returns:

       A generator object for Data frame(s) contains the browsing diseases information of metabolites from LMDB

    Raises:

        TypeError if argument (status) is not a list

        TypeError if argument (biospecimen) is not a list

        TypeError if argument (metabolite) is not a list

        TypeError if argument (disease) is not a list

    Example:

        for i_data in BroDis(status=["quantified"],biospecimen=["saliva","blood"],metabolite=["HMDB0002013","HMDB0002014"],disease=[],inborn_errors=True):
            print (i_data)

    """


    #make sure that status is fine to go
    if type(status) != list:
        raise TypeError ("status argument should be a list") 
    opp_status= ["quantified","detected","expected","predicted"]

    if len([i for i in status if i in opp_status]) != len(status):
        raise TypeError  ("status argument should be a list from [quantified, detected, expected ,or predicted]") 

    #make sure that biospecimen is fine to go
    if type(biospecimen) != list:
        raise TypeError ("biospecimen argument should be a list") 
    opp_biospecimen= ["blood","saliva","urine","csf","feces","sweat","breast_milk","bile","amniotic_fluid","other_fluids"]

    if len([i for i in biospecimen if i in opp_biospecimen]) != len(biospecimen):
        raise TypeError  ("biospecimen argument should be a list from [blood,saliva,urine,csf,feces,sweat,breast_milk,bile,amniotic_fluid,or other_fluids]") 

    #make sure that metabolite is fine to go
    if type(metabolite) != list:
        raise TypeError ("metabolite argument should be a list")

    #make sure that disease is fine to go
    if type(disease) != list:
        raise TypeError ("disease argument should be a list") 

    main_url = "http://www.hmdb.ca/diseases?"
    dict_url = {"utf8":"✓","filter":"true"}

    if status and len(status) != 0:
        for i_stat in status:
            dict_url[i_stat] = "1"
    if biospecimen and len(biospecimen) != 0:
        for i_bio in biospecimen:
            dict_url[i_bio] = "1"
    if metabolite and len(metabolite) !=0:
        dict_url["metabolite"]= "1"
        metabo_list = "\n".join(metabolite)
        dict_url["metabolite_list"] = metabo_list
    if disease and len(disease) !=0:
        dict_url["disease"]= "1"
        dis_list = "\n".join(disease)
        dict_url["disease_list"] = dis_list
    if inborn_errors and inborn_errors == True:
        dict_url["inborn_error"] = "1"
    

    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request 
    
    #get the diseases names form the xml
    xml = urllib.request.urlopen(api_request).read()
    
    def find_dis_name(xml):
        import re
        xml = str(xml)
        end_of_head = (hmdb.end() for hmdb in re.finditer('class="panel-heading"><strong>',xml))
        names = ( xml[i:i+100] for i in end_of_head )
        names_final = [name[:name.find("&nbsp;")] for name in names]
        
        return names_final

    try:
        tables = pd.read_html(api_request)
        diseases_names = find_dis_name(xml)

        for i_tables, name in zip(range(len(tables)), diseases_names):
            tables[i_tables].index.name = name

            tables_new = tables[i_tables]

            # to make the row.names =  to the id 
            index_name = range(0,len(tables_new))
            index_value = f"{diseases_names[0]}," *len(index_name)
            index_value = index_value.split(",")
            dictionary = dict(zip(index_name, index_value))
            tables_new = tables_new.rename(index= dictionary )

            yield tables_new

    except:
        yield "No diseases found"
        

def Biospecimen(status=list(),biospecimen=list(),metabolite=list(),disease=list()):
    """

    Browsing Bio specimen information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        status (list): list of keywords to filter by metabolite status (quantified, detected, expected, or/and predicted) (default all)

        biospecimen (list): list of keywords to filter by biospecimen (blood,saliva,urine,csf,feces,sweat,breast_milk,bile,amniotic_fluid, or/and other_fluids)

        metabolite_list (list): default = empty list. list of LMDB ID and/or metabolite names

        disease_list (list): default = empty list. list of disases names

    Returns:

       A generator object for Data frame(s) contains the browsing diseases information of metabolites from LMDB

    Raises:

        TypeError if argument (status) is not a list

        TypeError if argument (biospecimen) is not a list

        TypeError if argument (metabolite) is not a list

        TypeError if argument (disease) is not a list

    Example:

        for i_data in Biospecimen(status=["quantified"],biospecimen=["saliva","blood"],metabolite=["HMDB0002013","HMDB0002014"],disease=[]):
            print (i_data)

    """


    #make sure that status is fine to go
    if type(status) != list:
        raise TypeError ("status argument should be a list") 
    opp_status= ["quantified","detected","expected","predicted"]

    if len([i for i in status if i in opp_status]) != len(status):
        raise TypeError  ("status argument should be a list from [quantified, detected, expected ,or predicted]") 

    #make sure that biospecimen is fine to go
    if type(biospecimen) != list:
        raise TypeError ("biospecimen argument should be a list") 
    opp_biospecimen= ["blood","saliva","urine","csf","feces","sweat","breast_milk","bile","amniotic_fluid","other_fluids"]

    if len([i for i in biospecimen if i in opp_biospecimen]) != len(biospecimen):
        raise TypeError  ("biospecimen argument should be a list from [blood,saliva,urine,csf,feces,sweat,breast_milk,bile,amniotic_fluid,or other_fluids]") 

    #make sure that metabolite is fine to go
    if type(metabolite) != list:
        raise TypeError ("metabolite argument should be a list")

    #make sure that disease is fine to go
    if type(disease) != list:
        raise TypeError ("disease argument should be a list") 

    main_url = "http://www.hmdb.ca/biofluids?"
    dict_url = {"utf8":"✓","filter":"true"}

    if status and len(status) != 0:
        for i_stat in status:
            dict_url[i_stat] = "1"
    if biospecimen and len(biospecimen) != 0:
        for i_bio in biospecimen:
            dict_url[i_bio] = "1"
    if metabolite and len(metabolite) !=0:
        dict_url["metabolite"]= "1"
        metabo_list = "\n".join(metabolite)
        dict_url["metabolite_list"] = metabo_list
    if disease and len(disease) !=0:
        dict_url["disease"]= "1"
        dis_list = "\n".join(disease)
        dict_url["disease_list"] = dis_list

    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request 
    
    #get the diseases names form the xml
    xml = urllib.request.urlopen(api_request).read()
    
    def find_dis_name(xml):
        import re
        xml = str(xml)
        end_of_head = (hmdb.end() for hmdb in re.finditer('class="panel-heading"><strong>',xml))
        names = ( xml[i:i+100] for i in end_of_head )
        names_final = [name[:name.find("&nbsp;")] for name in names]
        
        return names_final

    try:
        tables = pd.read_html(api_request)
        diseases_names = find_dis_name(xml)

        for i_tables, name in zip(range(len(tables)), diseases_names):
            tables[i_tables].index.name = name

            tables_new = tables[i_tables]

            # to make the row.names =  to the id 
            index_name = range(0,len(tables_new))
            index_value = f"{diseases_names[0]}," *len(index_name)
            index_value = index_value.split(",")
            dictionary = dict(zip(index_name, index_value))
            tables_new = tables_new.rename(index= dictionary )

            yield tables_new

    except:
        yield "No diseases found"


def Bioclass(status=list(),biospecimen=list()):
    """

    Browsing diseases information about small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        status (list): list of keywords to filter by metabolite status (quantified, detected, expected, or/and predicted) (default all)

        biospecimen (list): list of keywords to filter by biospecimen (blood,saliva,urine,csf,feces,sweat,breast_milk,bile,amniotic_fluid, or/and other_fluids)

    Returns:

       A generator object for Data frame(s) contains the browsing diseases information of metabolites from LMDB

    Raises:

        TypeError if argument (status) is not a list

        TypeError if argument (biospecimen) is not a list

    Example:

        print (Bioclass(status=[],biospecimen=["other_fluids"]))

    """


    #make sure that status is fine to go
    if type(status) != list:
        raise TypeError ("status argument should be a list") 
    opp_status= ["quantified","detected","expected","predicted"]

    if len([i for i in status if i in opp_status]) != len(status):
        raise TypeError  ("status argument should be a list from [quantified, detected, expected ,or predicted]") 

    #make sure that biospecimen is fine to go
    if type(biospecimen) != list:
        raise TypeError ("biospecimen argument should be a list") 
    opp_biospecimen= ["blood","saliva","urine","csf","feces","sweat","breast_milk","bile","amniotic_fluid","other_fluids"]

    if len([i for i in biospecimen if i in opp_biospecimen]) != len(biospecimen):
        raise TypeError  ("biospecimen argument should be a list from [blood,saliva,urine,csf,feces,sweat,breast_milk,bile,amniotic_fluid,or other_fluids]") 


    main_url = "http://www.hmdb.ca/classyfication?"
    dict_url = {"utf8":"✓","filter":"true"}

    if status and len(status) != 0:
        for i_stat in status:
            dict_url[i_stat] = "1"
    if biospecimen and len(biospecimen) != 0:
        for i_bio in biospecimen:
            dict_url[i_bio] = "1"
    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request     
    try:
        all_tables = pd.DataFrame()
        page_count = 0
        while True:
            page_count += 1
            api_request = f"{api_request}&page={page_count}"
            tables = pd.read_html(api_request)
            all_tables = all_tables.append(tables[0], ignore_index=True, sort=False)
            if len(tables[0]) == 0:
                break

        all_tables.columns = all_tables.columns.droplevel(1)

        if len(all_tables) == 0:
            return "No diseases found"
        else:
            return all_tables.iloc[0:,:5]
    except:
        return "No diseases found"


def Search(query,searcher):
    """

    Searching for query on small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        query (str): a query string

        searcher (str): searcher to use from ( ["metabolites","diseases","pathways","proteins","reactions"] )

    Returns:

        list of HMDB IDs

    Raises:

        TypeError if argument (query) is not a string

        TypeError if argument (searcher) is not a string

    Example:
    
        data = Search("Dehydroepiandrosterone","metabolites")
        print (data)

    """
    
    if type(query) != str:
        raise TypeError ("query argument should be a string")
    if type(searcher) != str:
        raise TypeError ("searcher argument should be a string")
    opp_searcher = ["metabolites","diseases","pathways","proteins","reactions"]
    if searcher not in opp_searcher:
        raise TypeError (' searcher argument should be a string from ["metabolites","diseases","pathways","proteins","reactions"] ')

    main_url = "http://www.hmdb.ca/unearth/q?"
    dict_url = {"utf8":"✓","query":query,"searcher":searcher}
    
    api_request = urllib.parse.urlencode(dict_url)

    api_request = main_url + api_request

    xml = urllib.request.urlopen(api_request).read()

    def find_HMDB(xml):
        import re
        xml = str(xml)
        all_hmdb = (i.end() for i in re.finditer("HMDB",xml))
        all_hmdb = list( set ( [ f"HMDB{xml[i:i+7]}" for i in all_hmdb if xml[i:i+7].isdigit()  ] ) )
        return all_hmdb

    HMDB_ID = find_HMDB(xml)

    return HMDB_ID


def ChemQuery(start=100,end=200,search_type="molecular",filters=list()):
    """

    Search by molecular weight for small molecule metabolites found in the human body from The Human Metabolome Database (HMDB)

    Args:

        start (int): default = 100, the start of molecular weight to search from.

        end (int): default = 200, the end of molecular weight to search to.

        search_type: default = "molecular", either by the Molecular weight / Average mass (molecular) or by the Monoisotopic mass (monoisotopic)

        filters (list): list of keywords to filter with metabolite status (quantified, detected, or/and expected)

    Returns:

        Data frame contains the ChemQuery by molecular weight search result of metabolites from HMDB

    Raises:

        TypeError if argument (start) is not a int

        TypeError if argument (end) is not a int

        TypeError if argument (search_type) is not a string 

        TypeError if argument (filters) is not a list and out of this list [quantified, detected, or/and expected] 

    Example:

        for data in ChemQuery(start=100,end=105,search_type="molecular",filters=["quantified"]):
            print (data.head(8))
    
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
    main_url = "http://www.hmdb.ca/structures/search/metabolites/mass?"
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
        all_hmdb = (i.end() for i in re.finditer("HMDB",xml))
        all_hmdb = list( set ( [ f"HMDB{xml[i:i+7]}" for i in all_hmdb if xml[i:i+7].isdigit()  ] ) )
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

        Data frame object contains information retrieved from HMDB about metabolites returned from the Mass Spectrum search

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
    main_url= "http://www.hmdb.ca/spectra/ms/search?"
    dict_url= {"utf8":"✓","commit":"Search","database":"HMDB","ms_search_ion_mode":mode,"tolerance":tolerance,"tolerance_units":tolerance_unit}

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

    all_tables_all.sort_values(by=['Adduct','Monoisotopic Mass'],ascending=True,inplace=True)
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

        predicted : if "True" it will include predicted spectra. (default="False")

    Returns:

        Data frame object contains information retrieved from HMDB about metabolites returned from Tandom Mass Spectrum search

    Raises:

        TypeError if argument (masses) is not a list of length more than zero

    Example:

        result =LCMSMS(p_ion_mass=146.0,p_ion_tolerance=0.1,parent_ion_mass_tolerance_units="Da",ion_mode="positive",cid="low",peaks=["40.948 0.174","56.022 0.424"],mz_tolerance=0.5,mz_tolerance_units="Da",predicted=True)
        print (result)

   """

    main_url = "http://www.hmdb.ca/spectra/ms_ms/search?"
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
    elif predicted == "False":
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



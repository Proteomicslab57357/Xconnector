import urllib
import urllib.request
import pandas as pd 
from pandas import DataFrame

def BroComp(strain=list(),status=list(),compounds=list(),proteins=list(),pathways=list(),reactions=list()):

    """
    Browsing metabolites information about small molecule metabolites found in or produced by Saccharomyces cerevisiae (also known as Baker’s yeast and Brewer’s yeast) (YMDB).

    Args:

        strain (list): Filter by strain ["bakers_yeast","brewers_yeast"]

        status (list): Filter by status ["compound_quantified","compound_expected"]

        compounds (list): provide a list of YMDB IDs and/or metabolite names ("bakers_yeast","brewers_yeast")

        proteins (list): Provide a list of Uniprot IDs and/or protein names 

        pathways (list): Provide a list of SMPDB Pathwhiz IDs and/or pathway names 

        reactions (list): Provide a list of SMPDB Reaction IDs

    Returns:

        Data frame(s) contains general information retrieved from YMDB about metabolites returned from the search

    Raises:

        TypeError if argument (strain) is not a list

        TypeError if argument (status) is not a list

        TypeError if argument (compounds) is not a list
        
        TypeError if argument (proteins) is not a list

        TypeError if argument (pathways) is not a list

        TypeError if argument (reactions) is not a list


    Example:
    
        data = BroComp(strain=["bakers_yeast"],status=["compound_quantified","compound_expected"],compounds=["YMDB00001","YMDB00002"])
        print (data)

    """

    main_url = "http://www.ymdb.ca/compounds?"
    dict_url = {"utf8":"✓","filter":"true"}

    # make sure that strain is good to go
    if type(strain) != list:
        raise TypeError ("strain argument should be a list") 

    opp_strain = ["bakers_yeast","brewers_yeast"]
    if len( [ i for i in strain if i in opp_strain ] ) != len(strain):
        raise TypeError ("strain argument should be a list of ['bakers_yeast','brewers_yeast'] ") 


    #make sure that status is fine to go
    if type(status) != list:
        raise TypeError ("status argument should be a list")

    opp_status= ["compound_quantified","compound_expected"]
    if len([i for i in status if i in opp_status]) != len(status):
        raise TypeError  ("status argument should be a list from [quantified, detected, expected ,or predicted]") 

    #strain
    for i_strain in strain:
        dict_url[i_strain] = "1"

    for i_status in status:
        dict_url[i_status] = "1"
    
    #compunds

    if len(compounds) != 0:
        dict_url["compound"] = "1"
        dict_url["compound_list"] = "\n".join(compounds)

    #proteins

    if len(proteins) != 0:
        dict_url["protein"] = "1"
        dict_url["protein_list"] = "\n".join(proteins)
    
    #pathways

    if len(pathways) != 0:
        dict_url["pathway"] = "1"
        dict_url["pathway_list"] = "\n".join(pathways)

    #reaction

    if len(reactions) != 0:
        dict_url["reaction"] = "1"
        dict_url["reaction_list"] = "\n".join(reactions)


    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request

    tables = pd.read_html(api_request)
    tables = tables[0].drop(["Structure"],1)

    def formula(form):
        form = form[::-1]
        number = ""
        for i in range(len(form)):
            if form[i].isdigit() or form[i] == ".":
                number = number + form[i]
            else:
                break
        form_new = form.replace(number,"")
        
        return str( ( form_new[::-1], number[::-1] ) )

    for i_form in tables['Formula Weight']:
        i_form_new = formula(str(i_form))
        tables["Formula Weight"].replace([i_form], [i_form_new], inplace=True)

    return tables

def structures(accession):

    """

    Retrieve Structure image of small molecule metabolites found in or produced by Saccharomyces cerevisiae (also known as Baker’s yeast and Brewer’s yeast) (YMDB)


    Args:

        accessions (str): str of only one YMDB ID 

    Returns:

        png image for the metabolite structure

    Raises:

        TypeError if argument (accession) is not a string of only one YMDB ID

    Example:

        structures("YMDB00001")
   
    """

    if type(accession) != str:
        raise TypeError ("structures takes string as an argument with only one YMDB ID")

    api_request = f"http://www.ymdb.ca/structures/{accession}/image.png"
    image_name = f"{accession}.png"
    image= urllib.request.urlretrieve(api_request, image_name)

    return image

def Geninfo(accessions):
    """

    Retrieve information about small molecule metabolites found in or produced by Saccharomyces cerevisiae (also known as Baker’s yeast and Brewer’s yeast) (YMDB).

    Args:

        accessions (list): list of YMDB ID 

    Returns:

        A generator object for Data frame(s) contains general information retrieved from YMDB about metabolites

    Raises:

        TypeError if argument (accessions) is not a list

    Note:

        If an accession has no information in YMDB the function will return an empty data frame

    Example:

        for data in Geninfo(["YMDB00001","YMDB00002"]):
            print (data)

    """

    def add_all_df(dfdf_new):


        dfdf = dfdf_new.drop(["col1"], 1)
        dfdf = dfdf.transpose()
        dfdf.columns = list(dfdf_new["col1"])
        dfdf = dfdf.rename(index={"col2": dfdf["YMDB ID"][0]})

        return dfdf

    if type(accessions) != list:
        raise TypeError ("Geninfo takes list as an argument")

    main_url = "http://www.ymdb.ca/compounds/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        # get the data in form of tables
        try:
            try:
                tables = pd.read_html(api_request)
            except:
                tables =pd.read_html(api_request, header= 0 )
            #removing row with no data
            tables0 = tables[0].iloc[:,0:2]
            #rename the col 
            tables0.columns = ["col1", "col2"]
  
            #take only the general information
            take_from_tables0 = ["Version","Status","Creation Date","Update Date","YMDB ID","Common Name","Chemical Formula","Average Molecular Weight",
            "Monoisotopic Molecular Weight","IUPAC Name","Traditional Name","CAS Registry Number",
            "SMILES","InChI Identifier","InChI Key","Kingdom","Super Class","Class","Sub Class","Direct Parent",
            "Molecular Framework","Role","State","Cellular Locations","Biospecimen Locations","Tissue Locations",
            "DrugBank ID","FoodDB ID","Chemspider ID","KEGG Compound ID","ChEBI ID","PubChem Compound"]
            
            tables0 = tables0[tables0["col1"].isin(take_from_tables0)]
            #tables0.reset_index(drop = True, inplace=True)
            tables0 = add_all_df(tables0)
            yield tables0

        except:
            take_from_tables0 = ["Version","Status","Creation Date","Update Date","YMDB ID","Common Name","Chemical Formula","Average Molecular Weight",
            "Monoisotopic Molecular Weight","IUPAC Name","Traditional Name","CAS Registry Number",
            "SMILES","InChI Identifier","InChI Key","Kingdom","Super Class","Class","Sub Class","Direct Parent",
            "Molecular Framework","Role","State","Cellular Locations","Biospecimen Locations","Tissue Locations",
            "DrugBank ID","FoodDB ID","Chemspider ID","KEGG Compound ID","ChEBI ID","PubChem Compound"]
            
            not_found = ["NaN"] * len(take_from_tables0)
            not_found = dict( zip(take_from_tables0,not_found) )
            tables0 = pd.DataFrame(data=not_found , index= [0])

            tables0 = tables0.rename(index={0:i_acc })

            yield tables0


def ChemQuery(start=100,end=200,search_type="molecular"):
    """

    Search by molecular weight for small molecule metabolites found in or produced by Saccharomyces cerevisiae (also known as Baker’s yeast and Brewer’s yeast) (YMDB)

    Args:

        start (int): default = 100, the start of molecular weight to search from.

        end (int): default = 200, the end of molecular weight to search to.

        search_type: default = "molecular", either by the Molecular weight / Average mass (molecular) or by the Monoisotopic mass (monoisotopic)


    Returns:

        Data frame contains the ChemQuery by molecular weight search result of metabolites from YMDB

    Raises:

        TypeError if argument (start) is not a int

        TypeError if argument (end) is not a int

        TypeError if argument (search_type) is not a string 


    Example:

        for data in ChemQuery(start=100,end=105,search_type="molecular"):
            print (data.head(8))
    
    """

    #make sure that start and end is fine to go
    if type(start) != float or type(end) != float or start >= end:
        raise TypeError ("start and end arguments must be integer, and also the start must be less than the end")
    #make sure that search_type is molecular or monoisotopic
    if search_type != "molecular" and search_type != "monoisotopic":
        raise TabError ( "search_type argument must be either (molecular) for Molecular weight / Average mass, or (monoisotopic) for Monoisotopic mass" )


    #the main url
    main_url = "http://www.ymdb.ca/structures/search/compounds/mass?"
    dict_url = {"utf8":"✓","query_from":start,"query_to":end,"search_type":search_type,"commit":"Search"}

    #make the api url

    api_request =  urllib.parse.urlencode(dict_url)
    old_api_request = main_url + api_request


    def find_YMDB(xml):
        import re
        xml = str(xml)
        all_YMDB = (i.end() for i in re.finditer("YMDB",xml))
        all_YMDB = list( set ( [ f"YMDB{xml[i:i+5]}" for i in all_YMDB if xml[i:i+5].isdigit()  ] ) )
        return all_YMDB

    page = 0
    all_hmdb_id = []
    while True:
        api_request = old_api_request
        page = page + 1
        api_request = f"{api_request}&page={page}"
        xml = urllib.request.urlopen(api_request).read()
        final_YMDB = find_YMDB(xml)
        if len(final_YMDB) != 0:
            all_hmdb_id = all_hmdb_id + [ i for i in  final_YMDB] 
        if len(final_YMDB) == 0:
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

        Data frame object contains information retrieved from YMDB about metabolites returned from the Mass Spectrum search

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
    main_url= "http://www.ymdb.ca/spectra/ms/search?"
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

        predicted (str) : if True it will include predicted spectra. (default="False")

    Returns:

        Data frame object contains information retrieved from YMDB about metabolites returned from Tandom Mass Spectrum search

    Raises:

        TypeError if argument (masses) is not a list of length more than zero

    Example:

    result =LCMSMS(p_ion_mass=146.0,p_ion_tolerance=0.1,parent_ion_mass_tolerance_units="Da",ion_mode="positive",cid="low",peaks=["40.948 0.174","56.022 0.424"],mz_tolerance=0.5,mz_tolerance_units="Da",predicted=True)
    print (result)

   """

    main_url = "http://www.ymdb.ca/spectra/ms_ms/search?"
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


def txtsearch(query):
    """

    Searching using a powerful search engine based on the Lucene query language.

    Args:

        query (str): a query string

    Returns:

        A generator object for Data frame(s) contains information retrieved from LMDB about metabolites

    Raises:

        TypeError if argument (query) is not a string

    Example:
    
        data = txtsearch('"1-Methylhistidine"')
        print (data)

    """
    if type(query) != str:
        raise TypeError ("query argument should be a list")

    main_url = "http://www.ymdb.ca/unearth/q?"
    dict_url = {"utf8":"✓","query":query,"searcher":"compounds"}

    api_request = urllib.parse.urlencode(dict_url)
    api_request = main_url + api_request
    xml = urllib.request.urlopen(api_request).read()
  
    def find_YMDB(xml_str):
        import re
        xml_str = str(xml_str)
        if "returned no results" not in xml_str:
            end_of_lmdb = (lmdb.end() for lmdb in re.finditer("YMDB",xml_str))
            final_LMDB = list( set( [ f"YMDB{xml_str[end:end+5]}" for end in end_of_lmdb if xml_str[end:end+5].isdigit()] ) )
            return final_LMDB
        elif "returned no results" in xml_str:
            return False

    xml_parse =  find_YMDB(xml)

    return xml_parse

def ExpProp(accessions):

    """

    Retrieve Experimental Properties information about small molecule metabolites found in the human body from The Yeast Metabolome Database (YMDB)

    Args:

        accessions (list): list of YMDB ID

    Returns:

        A generator object for Data frame(s) contains the  Experimental Propertiess of metabolites from YMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in ExpProp(["YMDB00001","YMDB00002","YMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("ExpProp takes list as an argument")

    main_url = "http://www.ymdb.ca/compounds/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[1]
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

    Retrieve Predicted Properties information about small molecule metabolites found in the human body from The Yeast Metabolome Database (YMDB)

    Args:

        accessions (list): list of YMDB ID

    Returns:

        A generator object for Data frame(s) contains the Predicted Properties of metabolites from YMDB

    Raises:

        TypeError if argument (accessions) is not a list

    Example:

        for i_data in PredProp(["YMDB00001","YMDB00002","YMDB00003"]):
            print (i_data.head())

    """


    if type(accessions) != list:
        raise TypeError ("PredProp takes list as an argument")

    main_url = "http://www.ymdb.ca/compounds/"
    for i_acc in accessions:
        # creating the query search url by accessions
        api_request = main_url + i_acc
        try:
            tables = pd.read_html(api_request)
            tables = tables[2]
        except:
            tables = pd.DataFrame(data = {"Property":["NA"],"Value":["NA"],"Source":["NA"]})

        # to make the row.names =  to the id 
        index_name = range(0,len(tables))
        index_value = f"{i_acc}," *len(index_name)
        index_value = index_value.split(",")
        dictionary = dict(zip(index_name, index_value))
        tables = tables.rename(index= dictionary )

        yield tables

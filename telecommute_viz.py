##############
## Visualizing the future of remote work
##
## developed by Matt Delventhal
##
## Interactive visualization of results from "Spatial Implications of Telecommuting" (Delventhal and Parhomenko, 2022).
##
##
## For information on the research paper behind this app, please visit https://mattdelventhal.com/publication/spatial_implications_telecommuting/ 
##############

########
### Load dependencies

import math
import pandas as pd
import pydeck as pdk
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import streamlit as st

###
########


########
### Cached functions

## Sets color scale
@st.cache
def set_colorscale():
    breakz_0 = [-x for x in list(range(26))[::-1]] + list(range(26))[1:]
    curmapr = list(np.linspace(.9*.9,.975,1000))
    curmapr = curmapr[0:-1]
    curmapr += list(np.linspace(.975,.2*.8,1000))
    curmapg = list(np.linspace(.1*.8,.975,1000))
    curmapg = curmapg[0:-1];
    curmapg += list(np.linspace(.975,.1*.8,1000))
    curmapb = list(np.linspace(.1*.8,.975,1000))
    curmapb = curmapb[0:-1]
    curmapb += list(np.linspace(.975,.9*.8,1000))
    purplecoralmap3 = [[int(y*255) for y in x] for x in list(zip(curmapr,curmapg,curmapb))]

    colorrange_0 = [purplecoralmap3[int(x)] for x in np.linspace(0,len(purplecoralmap3)-1,len(breakz_0))]

    breakz_1 = [-x for x in list(range(11))[::-1]] + list(range(11))[1:]
    curmapr = list(np.linspace(.9*.9,.975,1000))
    curmapr = curmapr[0:-1]
    curmapr += list(np.linspace(.975,.1*.9,1000))
    curmapg = list(np.linspace(.3*.8,.975,1000))
    curmapg = curmapg[0:-1];
    curmapg += list(np.linspace(.975,.9*.9,1000))
    curmapb = list(np.linspace(.1*.9,.975,1000))
    curmapb = curmapb[0:-1]
    curmapb += list(np.linspace(.975,.4*.9,1000))
    greencoralmap = [[int(y*255) for y in x] for x in list(zip(curmapr,curmapg,curmapb))]

    colorrange_1 = [greencoralmap[int(x)] for x in np.linspace(0,len(greencoralmap)-1,len(breakz_1))]
    return [breakz_0,breakz_1],[colorrange_0,colorrange_1]

## Loads all the data
@st.cache(allow_output_mutation=True)
def load_process_main_data(breaklist,colorrangelist):
    
    # helper function to translate color range and breaks into a pydeck-compatible colorscale
    def color_scale(val,breakz,colorrange):
        for i, b in enumerate(breakz):
            if val < b:
                return colorrange[i]
        return colorrange[i]
    
    #####
    ## Load list of 50 largest U.S. cities with lat,lon to be selectable as map focii.
    df_state_to_select = pd.read_csv("states_to_select_withloc.csv")
    citylocdict = {}
    
    for j in df_state_to_select.index: #state in df_state_to_select["state"].unique():
        fullcityname = str(df_state_to_select["city"][j]) + ", " + str(df_state_to_select["state"][j])
        citylocdict[fullcityname] = {}
        citylocdict[fullcityname]["lat"] = float(df_state_to_select["lat"][j])
        citylocdict[fullcityname]["lon"] = float(df_state_to_select["lon"][j])
    ##
    #####
    
    ## Load model input and output data from Delventhal and Parkhomenko (2022).
    df_tele = pd.read_csv("result_data_export_aggl_altcf_increase_wfhutility.csv").rename(columns={"loccoder":"loccode"})
    
    ## Load county names, to match with county fips codes
    df_countynames = pd.read_csv("Countynames.csv")
    df_countynames['fullcounty'] = np.array(df_countynames.county) + np.array([', ']*len(df_countynames)) + np.array(df_countynames.state2dig)
    df_countynames['countycode'] = df_countynames['countycode'].astype('int64')

    ## Load PUMA names, to match with PUMA fips codes
    df_pumanames = pd.read_csv('2010_PUMA_Names.csv')
    
    ## Load shapefile containing polygons for model locations
    gdf_modellocs = gpd.read_file("model_locations_simplified_v10.shp")
    gdf_modellocs = gdf_modellocs.rename(columns={'lat':'lon','lon':'lat'})
    ### Simplify polygons to improve performance and make etruded shapes sleeker
    #gdf_modellocs["geometry"] = gdf_modellocs["geometry"].simplify(.004)

    #####
    ## Process shapefile polygons into a format acceptable to PyDeck.
    
    polygonlist = list(gdf_modellocs["geometry"])
    newpolylist = []
    loccodelist = []
    pumacodelist = []
    countycodelist = []
    latlist = []
    lonlist = []
    for j in range(len(polygonlist)):
        if type(polygonlist[j]) is Polygon:
            newpolylist.append(polygonlist[j])
            loccodelist.append(int(gdf_modellocs["loccode"][j]))
            pumacodelist.append(int(gdf_modellocs["pumacode"][j]))
            countycodelist.append(int(gdf_modellocs["countycode"][j]))
            latlist.append(float(gdf_modellocs["lat"][j]))
            lonlist.append(float(gdf_modellocs["lon"][j]))
        elif type(polygonlist[j]) is MultiPolygon:
            newpolylist += list(polygonlist[j])
            loccodelist += [int(gdf_modellocs["loccode"][j])]*len(list(polygonlist[j]))
            pumacodelist += [int(gdf_modellocs["pumacode"][j])]*len(list(polygonlist[j]))
            countycodelist += [int(gdf_modellocs["countycode"][j])]*len(list(polygonlist[j]))
            latlist += [float(gdf_modellocs["lat"][j])]*len(list(polygonlist[j]))
            lonlist += [float(gdf_modellocs["lon"][j])]*len(list(polygonlist[j]))
    
    allpolygons = [[list(x) for x in list(zip(poly.exterior.coords.xy[0],poly.exterior.coords.xy[1]))] for poly in newpolylist]
    ##
    #####
    
    #####
    ## Build dataframe with the processed polygons, and merge in the other data
    df_shapelayer = pd.DataFrame({'geometry':allpolygons,\
                                  'loccode':loccodelist,'pumacode':pumacodelist,'countycode':countycodelist,\
                                  'lat':latlist,'lon':lonlist})
    
    df_shapelayer = df_shapelayer.merge(df_tele[["loccode",
                                                 "sNRi_bl","sNRi_cf",
                                                 "sNWj_bl","sNWj_cf",
                                                 'Land','qhat']],how='left',left_on='loccode',right_on='loccode')
    df_shapelayer = df_shapelayer.merge(df_pumanames,how='left',left_on='pumacode',right_on='pumacode')
    df_shapelayer = df_shapelayer.merge(df_countynames,how='left',left_on='countycode',right_on='countycode')
    
    ##
    #####
    
    #####
    ## Check whether model location is built off of a county within a PUMA (rural), or a PUMA within a county (urban).
    ## For rural locations where there is nothing to say about the location other than what county it is, leave "n/a" for location detail.
    ## For urban location, provide location detail based on PUMA description.
    pumas_per_county = df_shapelayer.groupby(['countycode'])['pumacode'].nunique()
    df_shapelayer = df_shapelayer.merge(pumas_per_county.rename('pumas_per_county'),how='left',left_on='countycode',right_index=True)

    df_shapelayer['locdetail'] = 'n/a'
    df_shapelayer.loc[df_shapelayer.pumas_per_county > 1,'locdetail'] = df_shapelayer.loc[df_shapelayer.pumas_per_county > 1,'PUMANAME']
    ##
    #####
    
    #####
    ## Calculate basic variables: densities, absolute and relative changes, etc.
    df_shapelayer["sNRi_bl_dens"] = df_shapelayer["sNRi_bl"]/df_shapelayer['Land']
    df_shapelayer["sNRi_cf_dens"] = df_shapelayer["sNRi_cf"]/df_shapelayer['Land']
    df_shapelayer["Rihat_abs"] = df_shapelayer["sNRi_cf_dens"]  - df_shapelayer["sNRi_bl_dens"]
    df_shapelayer["Rihat_rel"] = df_shapelayer["Rihat_abs"]/df_shapelayer["sNRi_bl_dens"]

    df_shapelayer["sNWj_bl_dens"] = df_shapelayer["sNWj_bl"]/df_shapelayer['Land']
    df_shapelayer["sNWj_cf_dens"] = df_shapelayer["sNWj_cf"]/df_shapelayer['Land']
    df_shapelayer["Wjhat_abs"] = df_shapelayer["sNWj_cf_dens"] - df_shapelayer["sNWj_bl_dens"]
    df_shapelayer["Wjhat_rel"] = df_shapelayer["Wjhat_abs"]/df_shapelayer["sNWj_bl_dens"]
    ##
    ####
    
    #####
    ## Calculate variables to pass colors and heights to the map chart.
    df_shapelayer["Ri_elevation"] = df_shapelayer["sNRi_bl_dens"].apply(lambda x: (max(0,x-50))*1.5)
    df_shapelayer["Wj_elevation"] = df_shapelayer["sNWj_bl_dens"].apply(lambda x: (max(0,x-50))*1.5)
    df_shapelayer["q_elevation"] = df_shapelayer["sNRi_bl_dens"].apply(lambda x: (max(0,x-50))*1.5)

    df_shapelayer["Rihat_abs_c"] = df_shapelayer["Rihat_abs"].apply(lambda x: color_scale(posneg_sqrt(x),
                                                                                                        breaklist[0],
                                                                                                        colorrangelist[0]))
    df_shapelayer["Rihat_rel_c"] = df_shapelayer["Rihat_rel"].apply(lambda x: color_scale(100*x,
                                                                                                        breaklist[1],
                                                                                                        colorrangelist[1]))

    df_shapelayer["Wjhat_abs_c"] = df_shapelayer["Wjhat_abs"].apply(lambda x: color_scale(posneg_sqrt(x),
                                                                                                        breaklist[0],
                                                                                                        colorrangelist[0]))
    df_shapelayer["Wjhat_rel_c"] = df_shapelayer["Wjhat_rel"].apply(lambda x: color_scale(100*x,
                                                                                                        breaklist[1],
                                                                                                        colorrangelist[1]))
    
    df_shapelayer["qhat"] = df_shapelayer["qhat"] - 1
    df_shapelayer = df_shapelayer.rename(columns={'qhat':'qhat_rel'})
    df_shapelayer["qhat_rel_c"] = df_shapelayer["qhat_rel"].apply(lambda x: color_scale(100*x,
                                                                                                        breaklist[1],
                                                                                                        colorrangelist[1]))
    
    ##
    #####
    
    #####
    ## Calculate formatted variables, either for tooltip display, or for the user-downloadable data.    
    df_shapelayer["Ri_todisplay"] = df_shapelayer["sNRi_bl_dens"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["Ri_todisplayalt"] = df_shapelayer["sNRi_bl_dens"].apply(lambda x: '{:.1f}'.format(x))
    df_shapelayer["Rihat_abs_todisplay"] = df_shapelayer["Rihat_abs"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["Rihat_abs_todisplayalt"] = df_shapelayer["Rihat_abs"].apply(lambda x: '{:.1f}'.format(x))
    df_shapelayer["Rihat_rel_todisplay"] = df_shapelayer["Rihat_rel"].apply(lambda x: '{:.1f}'.format(100*x))
    df_shapelayer["Rihat_rel_todisplayalt"] = df_shapelayer["Rihat_rel"].apply(lambda x: '{:.4f}'.format(x))
    
    df_shapelayer["Wj_todisplay"] = df_shapelayer["sNWj_bl_dens"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["Wj_todisplayalt"] = df_shapelayer["sNWj_bl_dens"].apply(lambda x: '{:.1f}'.format(x))
    df_shapelayer["Wjhat_abs_todisplay"] = df_shapelayer["Wjhat_abs"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["Wjhat_abs_todisplayalt"] = df_shapelayer["Wjhat_abs"].apply(lambda x: '{:.1f}'.format(x))
    df_shapelayer["Wjhat_rel_todisplay"] = df_shapelayer["Wjhat_rel"].apply(lambda x: '{:.1f}'.format(100*x))
    df_shapelayer["Wjhat_rel_todisplayalt"] = df_shapelayer["Wjhat_rel"].apply(lambda x: '{:.4f}'.format(x))
        
    df_shapelayer["q_todisplay"] = df_shapelayer["sNRi_bl_dens"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["qhat_rel_todisplay"] = df_shapelayer["qhat_rel"].apply(lambda x: '{:.1f}'.format(100*x))
    df_shapelayer["qhat_rel_todisplayalt"] = df_shapelayer["qhat_rel"].apply(lambda x: '{:.4f}'.format(x))
    
    df_shapelayer["sNRi_bl_todisplay"] = df_shapelayer["sNRi_bl"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["sNWj_bl_todisplay"] = df_shapelayer["sNWj_bl"].apply(lambda x: '{:.0f}'.format(x))
    df_shapelayer["Land_todisplay"] = df_shapelayer["Land"].apply(lambda x: '{:.3f}'.format(x))
    ##
    #####
    
    #####
    ## Form the user-downloadable dataset. Give the variables nice, self-explanatory names.
    df_to_download = df_shapelayer[['loccode',
                                        'STATEFP',
                                        'state2dig',
                                        'countycode',
                                        'county',
                                        'pumacode',
                                        'PUMANAME',
                                        'sNRi_bl_todisplay',
                                        'sNWj_bl_todisplay',
                                        'Land_todisplay',
                                        'Ri_todisplayalt',
                                        'Wj_todisplayalt',
                                        'Rihat_abs_todisplayalt',
                                        'Rihat_rel_todisplayalt',
                                        'Wjhat_abs_todisplayalt',
                                        'Wjhat_rel_todisplayalt',
                                        'qhat_rel_todisplayalt']]
    
    df_to_download = df_to_download.rename(columns={'loccode':'unique_location_identifier',
                                        'STATEFP':'state_fips_code',
                                        'state2dig':'statename_abbrev',
                                        'countycode':'count_fips_code',
                                        'county':'countyname',
                                        'pumacode':'PUMA_fips_code',
                                        'PUMANAME':'PUMA_name',
                                        'sNRi_bl_todisplay':'residents_baseline',
                                        'sNWj_bl_todisplay':'jobs_baseline',
                                        'Land_todisplay':'landarea_sqkm',
                                        'Ri_todisplayalt':'residents_persqkm_baseline',
                                        'Wj_todisplayalt':'jobs_persqkm_baseline',
                                        'Rihat_abs_todisplay':'residents_persqkm_absolute_change',
                                        'Rihat_rel_todisplayalt':'residents_persqkm_relative_change',
                                        'Wjhat_abs_todisplay':'jobs_persqkm_absolute_change',
                                        'Wjhat_rel_todisplayalt':'jobs_persqkm_relative_change',
                                        'qhat_rel_todisplayalt':'floorspace_price_relative_change'})
    ##
    #####

    ## Drop unnecessary columns, on the off chance it may help performance.
    df_shapelayer.drop(columns=['Ri_todisplayalt','Wj_todisplayalt','Rihat_rel_todisplayalt','Wjhat_rel_todisplayalt',\
                                'Rihat_abs_todisplayalt','Wjhat_abs_todisplayalt',\
                                'qhat_rel_todisplayalt','sNRi_bl_todisplay','sNWj_bl_todisplay','Land_todisplay',\
                                'STATEFP','PUMA5CE','PUMANAME','pumas_per_county'])    
    
    ## Get rid of all "NaN"-coded Nan-values, to avoid a JSON error when the data gets passed to PyDeck.    
    for column in df_shapelayer.columns:
        if np.sum(np.array(df_shapelayer[column].isna())) > 0:
            if (df_shapelayer[column].dtype == "float64")\
            or (df_shapelayer[column].dtype == "float32")\
            or (df_shapelayer[column].dtype == "int")\
            or (df_shapelayer[column].dtype == "int8")\
            or (df_shapelayer[column].dtype == "int16")\
            or (df_shapelayer[column].dtype == "int32")\
            or (df_shapelayer[column].dtype == "int64"):
                df_shapelayer[column] = df_shapelayer[column].fillna(-9999)
            else:
                df_shapelayer[column] = df_shapelayer[column].fillna('')
    
    return df_shapelayer,citylocdict,df_to_download


## Converts pandas dataframe to a CSV ready to download.
@st.cache
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')

###
########


########
### Helper functions

## For x < 0, returns -sqrt(abs(x)).
## For x >= 0, returns sqrt(x).
def posneg_sqrt(x):
    if x < 0:
        return -math.sqrt(abs(x))
    else:
        return math.sqrt(x)

## Resets the map scrolling offsets. Used when "jumping" the view to a new city.
def reset_offset(targs):
    st.session_state[targs[0]] = 0
    st.session_state[targs[1]] = 0

## Define the dictionary of headings/bearings, needed for the setdirection function
bearingdict = {}
bearingdict['N'] = 0
bearingdict['NE'] = 45
bearingdict['E'] = 90
bearingdict['SE'] = 135
bearingdict['S'] = 180
bearingdict['SW'] = 225
bearingdict['W'] = 270
bearingdict['NW'] = 315

## Set the view direction in the session state.
def setdirection(direc,bearingdict=bearingdict):
    for key in list(bearingdict.keys()):
        if key == direc:
            st.session_state['direction'][key] = True
            st.session_state['curbearing'] = bearingdict[key]
        else:
            st.session_state['direction'][key] = False

###
########

########
### Set default values for scroll offset, view direction.

viewoffset = [0,0]
default_direc = 'N'

###
########


########
### Initialize session state variables
if 'direction' not in st.session_state:
    st.session_state['direction'] = {}
    st.session_state['curbearing'] = bearingdict[default_direc]
    for direc in list(bearingdict.keys()):
        if direc == default_direc:
            st.session_state['direction'][direc] = True
        else:
            st.session_state['direction'][direc] = False

if "viewoffset_lat" not in st.session_state:
    st.session_state["viewoffset_lat"] = 0

if "viewoffset_lon" not in st.session_state:
    st.session_state["viewoffset_lon"] = 0

###
########


########
### Define dictionaries to deploy options from radio buttons

vardict = {}
vardict['varselect'] = {}
vardict['varselect']["Change in residents"] = "Ri"
vardict['varselect']["Change in jobs"] = "Wj"
vardict['varselect']["Change in floorspace prices"] = "q"
vardict['fordisplay'] = {}
vardict['fordisplay']["Change in residents"] = "Change in residents"
vardict['fordisplay']["Change in jobs"] = "Change in jobs"
vardict['fordisplay']["Change in floorspace prices"] = "Change floorspace prices"
vardict['forheight'] = {}
vardict['forheight']["Change in residents"] = "Residents"
vardict['forheight']["Change in jobs"] = "Jobs"
vardict['forheight']["Change in floorspace prices"] = "Residents"
vardict['forchangeradio_i'] = {}
vardict['forchangeradio_i']["Change in residents"] = 0
vardict['forchangeradio_i']["Change in jobs"] = 0
vardict['forchangeradio_i']["Change in floorspace prices"] = 1
vardict['forchangeradio_d'] = {}
vardict['forchangeradio_d']["Change in residents"] = False
vardict['forchangeradio_d']["Change in jobs"] = False
vardict['forchangeradio_d']["Change in floorspace prices"] = True
vardict['forcheckbox'] = {}
vardict['forcheckbox']["Change in residents"] = "persons"
vardict['forcheckbox']["Change in jobs"] = "persons"
vardict['forcheckbox']["Change in floorspace prices"] = "residents"

changetypedict = {}
changetypedict['varselect'] = {}
changetypedict['varselect']["Absolute (persons/sq km)"] = "abs"
changetypedict['varselect']["Relative (%)"] = "rel"
changetypedict['fordisplay'] = {}
changetypedict['fordisplay']["Absolute (persons/sq km)"] = "/sq km"
changetypedict['fordisplay']["Relative (%)"] = " (%)"

###
########


########
### Call cache functions to calculate color scales and load and process all data

BREAKSlist,COLOR_RANGElist = set_colorscale() # Color scales

df_shapelayera,citylocdict,df_to_downloada = load_process_main_data(BREAKSlist,COLOR_RANGElist) # The main data work

fulldata_csv = convert_df(df_to_downloada) # Prepare tidy CSV file for download


###
########




########
### Deploy basic layout, later to be filled with content

st.write("### Predicted long-run changes due to increased remote work") # Heading
varselectcolumns = st.columns(3) # 3 columns for selecting city, variable to display, etc.

mainmap = st.container() # Container for main map

footer = st.container() # Footer

# Sidebar for map view tools
with st.sidebar:
    st.write("**Map view navigation**")
    compasscolumns_row1 = st.columns(3)
    compasscolumns_row2 = st.columns(3)
    compasscolumns_row3 = st.columns(3)

###
########





########
### Deploy components to upper 3 columns
focus_city = varselectcolumns[0].selectbox("Jump to city:",list(citylocdict.keys()),
                                            index=list(citylocdict.keys()).index("San Francisco, CA"),
                                            on_change=reset_offset,args=(["viewoffset_lat","viewoffset_lon"],))



var_choice = varselectcolumns[1].radio("Variable to color",["Change in residents","Change in jobs", "Change in floorspace prices"],index=0)

changetype_choice = varselectcolumns[2].radio("Type of change to color",["Absolute (persons/sq km)","Relative (%)"],\
                                                index=int(vardict['forchangeradio_i'][var_choice]),\
                                                disabled=bool(vardict['forchangeradio_d'][var_choice]))

###
########


########
### Deploy map view tools to sidebar

zoomit = compasscolumns_row2[1].slider("Zoom",3.0,12.0,value=8.5,step=0.5) # Zoom select tool, in the center of the tool panel (row,col) = (2,2)
scrollincrement = .4*(8.5/zoomit)**2 # Set the lat/lon increment for clicked scrolling according to zoom level.


## First row of buttons
if compasscolumns_row1[0].button('Scroll Northwest'): # Scroll to the northwest. Tool panel position (row, col) = (1,1).
    st.session_state["viewoffset_lat"] += .6*scrollincrement/math.sqrt(2) # Shortened scroll distance for diagonals, to avoid a disorienting jump.
    st.session_state["viewoffset_lon"] -= .6*scrollincrement/math.sqrt(2) 

compasscolumns_row1[0].checkbox('Face Northwest',value=st.session_state['direction']['NW'], #Face to the northwest. Tool panel position (row, col) = (1,1).
                               on_change=setdirection,args=('NW',))

if compasscolumns_row1[1].button('Scroll North'): # Scroll due north. Tool panel position (row, col) = (1,2).
    st.session_state["viewoffset_lat"] += scrollincrement

compasscolumns_row1[1].checkbox('Face North',value=st.session_state['direction']['N'], # Face due north. Tool panel position (row, col) = (1,2).
                               on_change=setdirection,args=('N',))

if compasscolumns_row1[2].button('Scroll Northeast'): # Scroll northeast. Tool panel position (row, col) = (1,3).
    st.session_state["viewoffset_lat"] += .6*scrollincrement/math.sqrt(2) # Shortened scroll distance for diagonals, to avoid a disorienting jump.
    st.session_state["viewoffset_lon"] += .6*scrollincrement/math.sqrt(2) 

compasscolumns_row1[2].checkbox('Face Northeast',value=st.session_state['direction']['NE'], #Face northeast. Tool panel position (row, col) = (1,2).
                               on_change=setdirection,args=('NE',))


## Second row of buttons
if compasscolumns_row2[0].button('Scroll West'): # Scroll due west. Tool panel position (row, col) = (2,1).
    st.session_state["viewoffset_lon"] -= scrollincrement

compasscolumns_row2[0].checkbox('Face West',value=st.session_state['direction']['W'], # Face due west. Tool panel position (row, col) = (2,1).
                               on_change=setdirection,args=('W',))


if compasscolumns_row2[2].button('Scroll East'): # Scroll due east. Tool panel position (row, col) = (2,3).
    st.session_state["viewoffset_lon"] += scrollincrement

compasscolumns_row2[2].checkbox('Face East',value=st.session_state['direction']['E'], # Face due east. Tool panel position (row, col) = (2,3).
                               on_change=setdirection,args=('E',))


## Second row of buttons
if compasscolumns_row3[0].button('Scroll Southwest'): # Scroll to the southwest. Tool panel position (row, col) = (3,1).
    st.session_state["viewoffset_lat"] -= .6*scrollincrement/math.sqrt(2) # Shortened scroll distance for diagonals, to avoid a disorienting jump.
    st.session_state["viewoffset_lon"] -= .6*scrollincrement/math.sqrt(2) 

compasscolumns_row3[0].checkbox('Face Southwest',value=st.session_state['direction']['SW'], # Face to the northwest. Tool panel position (row, col) = (3,1).
                               on_change=setdirection,args=('SW',))



if compasscolumns_row3[1].button('Scroll South'): # Scroll due south. Tool panel position (row, col) = (3,2).
    st.session_state["viewoffset_lat"] -= scrollincrement

compasscolumns_row3[1].checkbox('Face South',value=st.session_state['direction']['S'],
                               on_change=setdirection,args=('S',))

if compasscolumns_row3[2].button('Scroll Southeast'): # Scroll to the southeast. Tool panel position (row, col) = (3,3).
    st.session_state["viewoffset_lat"] -= .6*scrollincrement/math.sqrt(2) # Shortened scroll distance for diagonals, to avoid a disorienting jump.
    st.session_state["viewoffset_lon"] += .6*scrollincrement/math.sqrt(2) 

compasscolumns_row3[2].checkbox('Face Southeast',value=st.session_state['direction']['SE'], # Scroll to the southeast. Tool panel position (row, col) = (3,3).
                               on_change=setdirection,args=('SE',))

###
########


########
### Deploy content to footer


footer.download_button(
         label="Download underlying data + model predictions",
         data=fulldata_csv,
         file_name='DP2022_data_predictions.csv',
         mime='text/csv',
     )

footer.write(" ")

footer.write("""
 - *Baseline numbers are averages for 2012-2016.*
 - *Predictions reflect the March 2022 version of* [Spatial Implications of Telecommuting](https://mattdelventhal.com/publication/spatial_implications_telecommuting/ \"Delventhal and Parkhomenko (2022)\").
 - *For more details, please consult the research paper through the link above.*
"""
)

###
########



########
### Deploy the map figure

## Checkbox to select (or not) height variation to reflect density
extrusion=mainmap.checkbox("Show baseline " + vardict['forcheckbox'][var_choice] + "/sq km as height",value=True)

## Set viewstate based on user selections
view_state = pdk.ViewState(latitude=citylocdict[focus_city]["lat"] + st.session_state["viewoffset_lat"], \
       longitude= citylocdict[focus_city]["lon"] + st.session_state["viewoffset_lon"], \
       zoom= zoomit, minZoom=3,maxZoom= 12, pitch= 40, bearing= st.session_state['curbearing'],controller=True
)


## Build Pydeck polygon layer, selecting particular variables for coloring and/or height based on user selections.
lons = np.array(df_shapelayera['lon'])
lats = np.array(df_shapelayera['lat'])

selected = lats < citylocdict[focus_city]["lat"] + st.session_state["viewoffset_lat"] + 3.6*(10/zoomit)**2
selected &= lats > citylocdict[focus_city]["lat"] + st.session_state["viewoffset_lat"] - 3.6*(10/zoomit)**2
selected &= lons < citylocdict[focus_city]["lon"] + st.session_state["viewoffset_lon"] + 3.6*(10/zoomit)**2
selected &= lons > citylocdict[focus_city]["lon"] + st.session_state["viewoffset_lon"] - 3.6*(10/zoomit)**2

polygon_layer = pdk.Layer(
        "PolygonLayer",
        df_shapelayera.loc[selected,["geometry",
                        str(vardict['varselect'][var_choice]) + "_elevation",
                        str(vardict['varselect'][var_choice]) + "hat_" + str(changetypedict['varselect'][changetype_choice]) +"_c",
                        str(vardict['varselect'][var_choice]) + '_todisplay',
                        str(vardict['varselect'][var_choice]) + 'hat_' + str(changetypedict['varselect'][changetype_choice]) +'_todisplay',
                        'state2dig',
                        'county',
                        'locdetail']],
        id="df_shapelayera",
        opacity=0.7,
        stroked=True,
        get_polygon="geometry",
        filled=True,
        extruded=extrusion,
        wireframe=True,
        get_elevation= str(vardict['varselect'][var_choice]) + "_elevation",
        get_fill_color= str(vardict['varselect'][var_choice]) + "hat_" + str(changetypedict['varselect'][changetype_choice]) +"_c",
        get_line_color=[255, 255, 255],
        auto_highlight=True,
        pickable=True,
        highlightColor=[255,255,153,128],
    )

## Build subset of tooltip to report value of "height" variable
heightfield = ''
if True:  #extrusion: # Uncomment "extrusion:" to make display of "height" variable depend on the map actually charting that variable as polygon heights.
    heightfield = '<td><b>(height) ' + str(vardict['forheight'][var_choice]) + ' per sq km:</b></td>' +\
                   '<td>{' + str(vardict['varselect'][var_choice]) + '_todisplay}</td>' +\
                   '</tr>' +\
                   '<tr>' +\
                   '<td></td>' +\
                   '<td></td>' +\
                   '</tr>'

## Build main tooltip. This will by parsed by Pydeck, outputting html filled with variables from the Layer.
tooltip = {"html": '<table border = "1" style="width: 100%;font-size:10px">' +\
                   '<colgoup>' +\
                   '<col span="1" style="width: 60%;">' +\
                   '<col span="1" style="width: 40%;">' +\
                   '</colgroup>' +\
                   '<tbody>' +\
                   '<tr>' +\
                   '<td><b>(color) ' + str(vardict['fordisplay'][var_choice]) + str(changetypedict['fordisplay'][changetype_choice]) + ':</b></td>' +\
                   '<td>{' + str(vardict['varselect'][var_choice]) + 'hat_' + str(changetypedict['varselect'][changetype_choice]) +'_todisplay}</td>' +\
                   '</tr>' +\
                   '<tr>' +\
                    heightfield +\
                   '<tr>' +\
                   '<td><b>State:</b></td>' +\
                   '<td>{state2dig}</td>' +\
                   '</tr>' +\
                   '<tr>' +\
                   '<td><b>County:</b></td>' +\
                   '<td>{county}</td>' +\
                   '</tr>' +\
                   '<tr>' +\
                   '<td><b>Location detail:</b></td>' +\
                   '<td>{locdetail}</td>' +\
                   '</tr>' +\
                   '<tr>' +\
                   #'<td><b>Location code:</b></td>' +\
                   #'<td>{loccode}</td>' +\
                   #'</tr>' +\
                   '</tbody></table>',
          "style": {"width":"360px","font-family":"arial"}}

## Build the Deck. For some reason, Streamlit caching makes this slower.
r = pdk.Deck(
        layers=[polygon_layer],
        initial_view_state=view_state,
        map_style=pdk.map_styles.LIGHT,
        tooltip=tooltip
    )

## Deploy the Deck.
mainmap.pydeck_chart(r)


###
########




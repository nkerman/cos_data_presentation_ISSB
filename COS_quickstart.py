# %% [markdown]
# # A very intro guide to COS Data

# %% [markdown]
# ## Start by downloading some data
# * FUV
#   * Proposal 15366
#   * *Andy might recognize this proposal as his COS LP4 Spectral Resolution Program!*
# * NUV
#   * Picked at random

# %%
# Import the necessary libraries
from astroquery.mast import Observations # For actually searching and downloading
from pathlib import Path # Handling system paths
# For reading and editing fits astropy table format files
from astropy.table import Table
# For reading and editing general fits files
from astropy.io import fits
# For dealing with units and unit conversions
from astropy import units as u
# Plotting
import matplotlib.pyplot as plt
%matplotlib inline

# %% [markdown]
# ### We begin by downloading the final spectrum product (`X1DSUM`) of an FUV observation.

# %%
# Download an example FUV dataset using astroquery:
## Find all the observations from a single HST Proposal
obs_from_proposal = Observations.query_criteria(
    proposal_id="15366" # The Proposal ID of the observations to download
)
## Find all the data products for these observations
products_from_proposal = Observations.get_product_list(
    obs_from_proposal
)
print(f"Found {len(products_from_proposal)} total data products")
## Filter to a specific subset, i.e. of filetypes
products_to_download = Observations.filter_products(
    products_from_proposal,
    productSubGroupDescription=["X1DSUM","ASN"] # Filters to only the X1DSUM and ASN files
)
# Download the FUV products
download_table = Observations.download_products(products_to_download)
# Make a list of the local paths, aggregated by type of file
fuv_x1dsum_products = [Path(local_path) for local_path in download_table["Local Path"] if "x1dsum" in local_path]
fuv_asn_products = [Path(local_path) for local_path in download_table["Local Path"] if "asn" in local_path]
print("FUV X1DSUM Files: \n", fuv_x1dsum_products, "\nFUV ASN Files: \n", fuv_asn_products)

# %% [markdown]
# ### We also want to download an example of the raw data (`RAWTAG`), and the intermediate data step (`CORRTAG`)
# This is not all the data which went into the final products downloaded above.
# It corresponds to a single exposure on a single segment of the FUV detector.

# %%
rawtag_a = Path(Observations.download_products(
    Observations.filter_products(
        Observations.get_product_list(
            Observations.query_criteria(
                proposal_id="15366"
            )
        ),
        productSubGroupDescription=["RAWTAG_A"]
    )[0]
)["Local Path"][0])

corrtag_a = Path(Observations.download_products(
    Observations.filter_products(
        Observations.get_product_list(
            Observations.query_criteria(
                proposal_id="15366"
            )
        ),
        productSubGroupDescription=["CORRTAG_A"]
    )[0]
)["Local Path"][0])

print(f"\n\nRaw TIME-TAG data from segment A in: {rawtag_a}")
print(f"Corrected TIME-TAG data from segment A in: {corrtag_a}")

# %% [markdown]
# ### Finally, let's download an example NUV dataset using astroquery in a more condensed form

# %%
nuv_x1dsum_g230l = Path(Observations.download_products(
    Observations.filter_products(
        Observations.get_product_list(
            Observations.query_criteria(
                obs_id = 'lbbd01020'
            )
        ),
        productSubGroupDescription=["X1DSUM"]
    )[0]
)["Local Path"][0])

nuv_x1dsum_g185m = Path(Observations.download_products(
    Observations.filter_products(
        Observations.get_product_list(
            Observations.query_criteria(
                obs_id = 'LAAD020L0'
            )
        ),
        productSubGroupDescription=["X1DSUM"]
    )[0]
)["Local Path"][0])

# %% [markdown]
# ## Examining the data products
# 
# ### The RAW and CORR Data:

# %%
hdr0 = fits.getheader(rawtag_a)
hdr1 = fits.getheader(rawtag_a, ext=1)
print(f"This is a {hdr1['EXPTIME']} second exposure on segment {hdr0['SEGMENT']} with the {hdr0['OPT_ELEM']} grating and the {hdr0['CENTRWV']} Å central wavelength setting.")

# My favorite way to read in FITS data - there are others!
raw_data = Table.read(rawtag_a, hdu=1)
corr_data = Table.read(corrtag_a, hdu=1)

# %%
raw_data

# %%
corr_data

# %% [markdown]
# ### Let's take a look at how the counts fell on the detector.

# %%
# So this runs quickly, limit to first 50k counts

plt.figure(figsize=(12,8))
plt.scatter(raw_data["RAWX"][:50000], raw_data["RAWY"][:50000], s=1, alpha=0.5, label="RAW Data")
plt.scatter(corr_data["XFULL"][:50000], corr_data["YFULL"][:50000], s=1, alpha=0.8, label="Corrected Data")
plt.xlabel("Detector X Coordinates")
plt.ylabel("Detector Y Coordinates")
plt.legend()
plt.show()

# %% [markdown]
# The power of this type of "`TIME-TAG`"data is that we have more control and information on each individual photon encounter. Because of this we have the ability to...
# * Remove suspected noise
# * Filter out "bad" time intervals
#   * i.e. to remove data taken when the sun was "up" for HST [see DayNight.ipynb](https://spacetelescope.github.io/COS-Notebooks/DayNight.html)
# * Split apart our data and examine it at different times
#   * i.e. for a transit [see SplitTag.ipynb](https://spacetelescope.github.io/COS-Notebooks/SplitTag.html)

# %% [markdown]
# ### Now let's examine the FUV spectrum

# %%
for fuv_x1d in fuv_x1dsum_products:
    hdr0 = fits.getheader(fuv_x1d)
    hdr1 = fits.getheader(fuv_x1d, ext=1)
    print(f"{hdr0['ASN_ID']} is a {hdr1['EXPTIME']} second (combined) exposure with the {hdr0['OPT_ELEM']} grating and the {hdr0['CENTRWV']} Å central wavelength setting.")

# %% [markdown]
# ### We'll look at the structure of the first (G130M) spectrum:
# * We note there are two rows for each of the segments of the detector
#   * If we want to plot, we'll need to deal with each of these segments separately
# * For NUV data, there are three rows for each "Data stripe"

# %%
fuv_g130m = Table.read(fuv_x1dsum_products[0])
fuv_g130m

# %%
# Set up the plot as a single box with size of 10x4 inches, and with a dpi of 100, relevant should we choose to save it:
fig1, ax = plt.subplots(1,1,figsize=(10,4), dpi = 100)  

# The next few lines are the core of the cell, where we actually place the data onto the plot:
###############
# Access each row of data, the longer wvln segment and gets the data we need to plot a spectrum:
for i, row in enumerate(fuv_g130m):
    wvln, flux, segment  = row["WAVELENGTH", "FLUX", "SEGMENT"] 
    ax.plot(wvln, flux, # First two arguments are assumed to be the x-data, y-data
            linestyle = "-", linewidth = 0.25, c = "km"[i], # These parameters specify the look of the connecting line
            marker = '.', markersize = 2, markerfacecolor = 'r', markeredgewidth = 0, # The marker parameters specify how the data points will look... 
                                                                                    # ... if you don't want dots set marker = ''
            label = segment) # The label is an optional parameter which will allow us to create a legend 
                            # this label is useful when there are multiple datasets on the same plot
# The lines after this are all about formatting, adding text, and saving as an image
###############
ax.set_title("G130M COS Spectrum", size = 20) # Adds a title of fontsize 20 points
ax.set_xlabel('Wavelength [$\AA$]', size = 12) # Adds x axis label
ax.set_ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size = 12) # Adds y label

plt.legend(loc = 'upper right') # Adds a legend with the label specified in the plotting call
plt.tight_layout() # Trims blank space
plt.show() # Shows all the plot calls in this cell and "clears" the plotting space - must come after any saving you want to do

# %%
# Set up the plot as a single box with size of 10x4 inches, and with a dpi of 100, relevant should we choose to save it:
fig1, ax = plt.subplots(1,1,figsize=(10,4), dpi = 100)  

# The next few lines are the core of the cell, where we actually place the data onto the plot:
###############
for j, fuv_file in enumerate(fuv_x1dsum_products):
    fuv_tab = Table.read(fuv_file)
    data_label = f"{fits.getval(fuv_file, 'OPT_ELEM')}/{fits.getval(fuv_file, 'CENTRWV'):.0f} "
    for i, row in enumerate(fuv_tab):
        wvln, flux, segment  = row["WAVELENGTH", "FLUX", "SEGMENT"] 
        ax.plot(wvln, flux, # First two arguments are assumed to be the x-data, y-data
                linestyle = "-", linewidth = 0.5, c = ["bg","mk","yc"][j][i], # These parameters specify the look of the connecting line
                marker = '', markersize = 2, markeredgewidth = 0, # The marker parameters specify how the data points will look... 
                label = data_label+segment, alpha=1-0.1*j)
###############
ax.set_title("G130M+G160M COS Spectrum", size = 20) # Adds a title of fontsize 20 points
ax.set_xlabel('Wavelength [$\AA$]', size = 12) # Adds x axis label
ax.set_ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size = 12) # Adds y label

plt.legend(loc = 'upper right') # Adds a legend with the label specified in the plotting call
plt.tight_layout() # Trims blank space
plt.show() # Shows all the plot calls in this cell and "clears" the plotting space - must come after any saving you want to do

# %% [markdown]
# ### Finally let's look at some NUV spectra
# Both spectra are of white dwarf stars
# 1. First the G185M Grating spectrum of GD71
# 2. Second the G230L Grating spectrum of WD1057+719
#    1. This is a somewhat strange spectrum COS creates with the G230L grating. The first two stripes are (like all other COS spectra) the first order spectrum of the source. The third stripe is dominated by second order light. It's thus more dispersed, has lower S/N, and falls on top of the first (NUVA) stripe.

# %%
nuv_tab = Table.read(nuv_x1dsum_g185m)
# Set up the plot as a single box with size of 10x4 inches, and with a dpi of 100, relevant should we choose to save it:
fig1, ax = plt.subplots(1,1,figsize=(10,4), dpi = 100)  

# The next few lines are the core of the cell, where we actually place the data onto the plot:
###############
# Access each row of data, the longer wvln segment and gets the data we need to plot a spectrum:
for i, row in enumerate(nuv_tab):
    wvln, flux, segment  = row["WAVELENGTH", "FLUX", "SEGMENT"] 
    ax.plot(wvln, flux, # First two arguments are assumed to be the x-data, y-data
            linestyle = "-", linewidth = 0.25, c = "kmb"[i], # These parameters specify the look of the connecting line
            marker = '.', markersize = 2, markerfacecolor = 'r', markeredgewidth = 0, # The marker parameters specify how the data points will look... 
            label = segment)
###############
ax.set_title("G185M NUV COS Spectrum", size = 20) # Adds a title of fontsize 20 points
ax.set_xlabel('Wavelength [$\AA$]', size = 12) # Adds x axis label
ax.set_ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size = 12) # Adds y label

plt.legend(loc = 'upper right') # Adds a legend with the label specified in the plotting call
plt.tight_layout() # Trims blank space
plt.show() # Shows all the plot calls in this cell and "clears" the plotting space - must come after any saving you want to do

# %%
nuv_tab = Table.read(nuv_x1dsum_g230l)
# Set up the plot as a single box with size of 10x4 inches, and with a dpi of 100, relevant should we choose to save it:
fig1, ax = plt.subplots(1,1,figsize=(10,4), dpi = 100)  

# The next few lines are the core of the cell, where we actually place the data onto the plot:
###############
# Access each row of data, the longer wvln segment and gets the data we need to plot a spectrum:
for i, row in enumerate(nuv_tab):
    wvln, flux, segment  = row["WAVELENGTH", "FLUX", "SEGMENT"] 
    ax.plot(wvln, flux, # First two arguments are assumed to be the x-data, y-data
            linestyle = "-", linewidth = 0.25, c = "kmb"[i], # These parameters specify the look of the connecting line
            marker = '.', markersize = 2, markerfacecolor = (1,0,0,0.5), markeredgewidth = 0, # The marker parameters specify how the data points will look... 
            label = segment)
###############
ax.set_title("G230L NUV COS Spectrum", size = 20) # Adds a title of fontsize 20 points
ax.set_xlabel('Wavelength [$\AA$]', size = 12) # Adds x axis label
ax.set_ylabel('Flux [$erg\ s^{-1}\ cm^{-2}\ Angstrom^{-1}$]', size = 12) # Adds y label

plt.legend(loc = 'upper right') # Adds a legend with the label specified in the plotting call
plt.tight_layout() # Trims blank space
plt.show() # Shows all the plot calls in this cell and "clears" the plotting space - must come after any saving you want to do

# %%




import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import matplotlib as mpl
import ftplib


#Todo: WIE GEHT DAS!?!??!?!?!?!?!?!!?
#meinftp = ftplib.FTP("ftp://arthurhou.pps.eosdis.nasa.gov")

#meinftp.login("bregovic@gmx.de","bregovic@gmx.de")



A = "3B-HHR.MS.MRG.3IMERG.20150404-S000000-E002959.0000.V03D.HDF5"

ftp = ftplib.FTP("ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/2015/04/04/imerg/")
ftp.login("bregovic@gmx.de", "bregovic@gmx.de")
ftp.cwd("/Dir")





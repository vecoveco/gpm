"""
-----------------------------------------------
Program zum Downloaden von GPM Daten vom Server
-----------------------------------------------
"""

from pcc import melde_dich
import ftplib
from datetime import date, timedelta as td
import pandas as pd


meinftp = ftplib.FTP("arthurhou.pps.eosdis.nasa.gov")

meinftp.login("bregovic@gmx.de","bregovic@gmx.de")


a = pd.read_csv('/automount/ags/velibor/text/overpass_dpr_radolan.csv', sep=',')
overpass_zeiten = a[[0]].values

#zum auslasten 0,500; 500,1000; 1000,1500

for i in range(len(overpass_zeiten)-1):

    overpass_datum = overpass_zeiten[i][0][0:8]

    dpr_overpass = '2A.GPM.DPR.V7-20170308.' +overpass_zeiten[i][0][:-7]+'5A.HDF5'


    print overpass_zeiten[i][0][:-7]+'5A.HDF5'


    print dpr_overpass

    # DPR gibt es erst ab 2014/03/09

    directory = '/gpmdata/'+ str(overpass_datum[0:4]) + '/' + str(overpass_datum[4:6]) + '/' + str(overpass_datum[6:8]) +'/radar/'

    meinftp.cwd(directory)

    # Mein Verzeichniss
    directory_local = '/automount/ags/velibor/gpmdata/dprV7/'


    #meinftp.retrlines('LIST')

    daten=[]           #Initialisierung einer Liste (leere Liste)

    meinftp.dir(daten.append)


    print
    print 'Ort und Name der lokalen Datei: ' + directory_local + dpr_overpass
    print

    file = open(directory_local+dpr_overpass, 'wb')

    print 'Download: ftp-Server: ' + directory +dpr_overpass

    meinftp.retrbinary('RETR '+dpr_overpass, file.write)

    print
    print 'Die lokale Datei ' + directory_local+dpr_overpass +' wird geschlossen.'


    file.close()

print meinftp.quit()
print
print 'Die FTP-Verbindung wurde von mir getrennt.'


melde_dich('Das Program ist fertig!')
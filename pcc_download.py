"""
-----------------------------------------------
Program zum Downloaden von GPM Daten vom Server
-----------------------------------------------
"""


import ftplib
from datetime import date, timedelta as td

meinftp = ftplib.FTP("arthurhou.pps.eosdis.nasa.gov")

meinftp.login("bregovic@gmx.de","bregovic@gmx.de")



d1 = date(2015, 1, 1)
d2 = date(2015, 12, 31)

delta = d2 - d1

for i in range(delta.days + 1):
    zeit = d1 + td(days=i)
    print zeit

    directory = '/gpmdata/'+ str(zeit.strftime("%Y/%m/%d")) +'/gprof/'

    meinftp.cwd(directory)

    # Mein Verzeichniss
    directory_local = '/automount/ags/velibor/gpmdata/gprof/'


    #meinftp.retrlines('LIST')

    daten=[]           #Initialisierung einer Liste (leere Liste)

    meinftp.dir(daten.append)

    gprof_wort = '2A.GPM.GMI.GPROF2014v2-0.'

    for i in range(len(daten)):

        if gprof_wort in daten[i][-66::]:

            print '>>>>>>>> ', daten[i][-66::]

            filename1 = daten[i][-66::]
            filename2 = daten[i][-66::]

            print
            print 'Ort und Name der lokalen Datei: ' + directory_local + filename2
            print

            file = open(directory_local+filename2, 'wb')

            print 'Download: ftp-Server: ' + directory +filename1

            meinftp.retrbinary('RETR '+filename1, file.write)

            print
            print 'Die lokale Datei ' + directory_local+filename2 +' wird geschlossen.'


            file.close()

print meinftp.quit()
print
print 'Die FTP-Verbindung wurde von mir getrennt.'


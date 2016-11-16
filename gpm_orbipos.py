import ftplib
import urllib


my_host = 'ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/geolocation/2015/04/10/'
my_pass = 'bregovic@gmx.de'
my_usr = 'bregovic@gmx.de'



urllib.urlretrieve('ftp://bregovic@gmx.de:bregovic@gmx.de@ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/geolocation/2015/04/10/',
                   'GPMCORE.20150410.005816815_20150410.023048999.001.EPHEM.txt')
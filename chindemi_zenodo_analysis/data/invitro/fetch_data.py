#!/usr/bin/env python
# coding: utf-8
"""
Script to download and reformat the distance dependent data of 
Sjostrom and Hausser 2006

Running this script will generate the file 
"epsp_vs_dist_sjostrom_hausser_2006.csv" in the current working directory

Requirements/Dependencies:
$ pip install xlrd wget pandas
"""
import pandas as pd
import wget
import hashlib

data_url = "http://plasticity.muhc.mcgill.ca/DataPage/Sjostrom_2006_Fig_3/LTP_vs_loc.xls"
wget.download(data_url)

#Check the file integrity
check_file = hashlib.md5(open('LTP_vs_loc.xls','rb').read()).hexdigest()
if not check_file=="98e3085f0afe217acfff62a6cd4f9a33":
    raise ValueError("The downloaded data appears to be corrupt or altered.")

print("\n")

data = pd.read_excel (r'LTP_vs_loc.xls') 
data["type"] = "none"

l5 = pd.DataFrame(data, columns= ["type","Plast","RT"]).dropna()
l5["type"]="l5"

l23 = pd.DataFrame(data, columns= ["type","Plast.1","RT.1"]).dropna()
l23 = l23.rename(columns={"Plast.1":"Plast", "RT.1":"RT"})
l23["type"] = "l23"

ext = pd.DataFrame(data, columns= ["type","Plast.2","RT.2"]).dropna()
ext = ext.rename(columns={"Plast.2":"Plast", "RT.2":"RT"})
ext["type"] = "ext"

df = pd.concat([l5,l23,ext])
df.to_csv("plast_vs_dist_sjostrom_hausser_2006.csv",index=False, header=["stim", "epsp_ratio","risetime"])
print("Success.")





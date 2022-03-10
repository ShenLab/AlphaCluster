from os import walk
import os
import mygene
import sys


transcript_list_file=sys.argv[1]
with open(transcript_list_file) as f:
    filenames = f.readlines()
    #filenames = next(walk("../coevol"), (None, None, []))[2]
#transcript_name = sys.argv[1]
mg = mygene.MyGeneInfo()

for transcript_name in filenames:
    print(transcript_name)
    if 'ENST' in transcript_name:   
        #transcript_id = transcript_name.split('.')[2].split('/')[2]
        transcript_id = transcript_name.split('.')[0]
        print(transcript_id)
        symbol = mg.query("ensembl.transcript:"+transcript_id, field=["symbol"])
        if not len(symbol["hits"]) == 0 and "symbol" in symbol["hits"][0]:
            print(symbol["hits"][0]["symbol"])
            new_name = symbol["hits"][0]["symbol"]
            if not new_name == "":
                os.rename("../coevol/" + transcript_id+".coevol.txt", "../coevol/" + new_name +".coevol.txt")

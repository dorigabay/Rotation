import re
import os
import glob


def edit_MDwrapper(directory,MDwrapper_origin_file,output_path,segment,temp,salt):
    MDwrapper = open(MDwrapper_origin_file,"r")
    raw_text = MDwrapper.read().split("\n")
    new_MDwrapper = open(directory+"MDWrapper.prefs","w")
    objects_to_modify = {
        "OUTPUT_PATH":"",
        "IONIC_STRENGTH":"",
        "INPUT_FILE_PATH":"",
        "segA":"1,151",
        "segB":"1,159",
        "segC":"1,128",
        "TEMPERATURES":""
    }
    objects_to_modify["IONIC_STRENGTH"] = salt
    objects_to_modify["TEMPERATURES"] = temp
    objects_to_modify["OUTPUT_PATH"]=output_path
    print(objects_to_modify["OUTPUT_PATH"])
    objects_to_modify["INPUT_FILE_PATH"]= "new_file.dat"
    new_MDwrapper.write(raw_text[0])
    for row in raw_text[1:]:
        if re.match("OUTPUT_PATH.*",row):
            new_MDwrapper.write("\n"+"OUTPUT_PATH=/cluster_data/dor/coarse_grain/"+objects_to_modify["OUTPUT_PATH"])
        elif re.match("IONIC_STRENGTH.*", row):
            new_MDwrapper.write("\n"+"IONIC_STRENGTH="+objects_to_modify["IONIC_STRENGTH"])
        elif re.match("INPUT_FILE_PATH.*", row):
            new_MDwrapper.write("\n"+"INPUT_FILE_PATH="+"./"+objects_to_modify["INPUT_FILE_PATH"])
        elif re.match("DYNAMIC_ATOMS_RANGE.*", row):
            new_MDwrapper.write("\n"+"DYNAMIC_ATOMS_RANGE="+objects_to_modify[segment])
        elif re.match("TEMPERATURES.*", row):
            new_MDwrapper.write("\n"+"TEMPERATURES="+objects_to_modify["TEMPERATURES"])
        else:
            new_MDwrapper.write("\n"+row)


def iter_edit_MDwrapper(directory):
    directories_toiter = [directory + "/"]
    origin_MDwrapper = directory + "/MDWrapper.prefs"
    while directories_toiter:
        dir = directories_toiter.pop()
        if glob.glob(dir + "/*.dat"):
            segment = re.search(".*(seg[ABC]).*", dir).group(1)
            salt = str(int(re.search(".*salt(.*)mM.*",dir).group(1))/1000)
            print(salt)
            output_inolympus = re.search(directory+"(.*)",dir).group(1)
            edit_MDwrapper(dir, origin_MDwrapper, output_inolympus, segment, "0.4", salt)
        else:
            for subdir in [x for x in os.walk(dir)][0][1]:
                directories_toiter.append(dir + subdir + "/")


def create_empty_dirs(directory):
    directories_toiter = [directory + "/"]
    while directories_toiter:
        dir = directories_toiter.pop()
        if glob.glob(dir + "/*.dat"):
            os.makedirs(directory+"/directory_to_cluster/"+re.search(directory+"(.*)",dir).group(1))
        else:
            for subdir in [x for x in os.walk(dir)][0][1]:
                directories_toiter.append(dir + subdir + "/")

# iter_edit_MDwrapper("/home_e/dor/CPEB4/coarse_grain/HIS_charges")
create_empty_dirs("/home_e/dor/CPEB4/coarse_grain/HIS_charges")
import glob
import re
import os

seqs = {"seqA": "GSGGGSGGGSG",
        "seqB": "GSGGPSGGGSG",
        "seqC": "GSGPPSGGGSG",
        "seqD": "GSGGPSPGGSG",
        "seqE": "GSGGPSGPGSG",
        "seqF": "GSGGPPPGGSG",
        "seqG": "GSGPSPGPGSG",
        "seqH": "GSGGPPSPGSG"}

segment_seqs = {"segA":"CNTGNKSAFPVRFHPHLQPPHHHQNATPSPAAFINNNTAANGSSAGSAWLFPAPATHNIQDEILGSEKAKSQQQEQQDPLEKQQLSPSPGQEAGILPETEKAKSEENQGDNSSENGNGKEKIRIESPVLTGFDYQEATGLGTSTQPLTSSC",
                "segB":"CSSLTGFSNWSAAIAPSSSTIINEDASFFHQGGVPAASANNGALLFQNFPHHVSPGFGGSFSPQIGPLSQHHPHHPHFQHHHSQHQQQRRSPASPHPPPFTHRNAAFNQLPHLANNLNKPPSPWSSYQSPSPTPSSSWSPGGGGYGGWGGSQGRDHRRC",
                "segC":"CLNGGITPLNSISPLKKNFASNHIQLQKYARPSSAFAPKSWMEDSLNRADNIFPFPDRPRTFDMHSLESSLIDIMRAENDTIKARTYGRRRGQSSLFPMEDGFLDDGRGDQPLHSGLGSPHSFSHQNC"}



class __Drule1__:
    combi = ["([^P]{2}P.)", "P[^P](P[^P]P[^P])"]
    angel1 = "4.040"
    strength1 = "1.000"
    angel2 = "0.000"
    strength2 = "0.000"
#
#
class __Drule2__:
    combi = ["[^P].(P[^P]P.)", "([^P]PP.)", "(PPP.)", "P([^P]P[^P]P)","(PP[^P]P)"]
    angel1 = "4.890"
    strength1 = "2.000"
    angel2 = "0.000"
    strength2 = "0.000"
#
#
class __Drule3__:
    combi = ["PP([^P]P[^P][^P])[^P]", "(?=[^P][^P]P([^P]P[^P][^P]))", "P(PP[^P][^P])[^P]", "P[^P]P([^P]P[^P][^P])[^P]",
             "P[^P]([^P]P[^P][^P])","[^P][^P]([^P]P[^P][^P])[^P][^P]"]
    angel1 = "5.140"
    strength1 = "1.500"
    angel2 = "0.000"
    strength2 = "0.000"
#
#
class __Prule1__:
    combi = ["(?=(.P.))","(?=(PP.))","(?=(.PP))","(?=(P.P))"]
    angel1 = "1.200"
    strength1 = "10.000"

class __Prule2__:
    combi = ["(?=([^P]{3}))","(?=(P[^P][^P]))","(?=([^P][^P]P))"]
    angel1 = "0.800"
    strength1 = "10.000"


def P_search_rule(pro_seq):
    angel_dict = {__Prule1__:[],__Prule2__:[]}
    for combi in __Prule1__.combi:
        it = re.finditer(combi, pro_seq)
        locations = [x.start() for x in it]
        for i in locations:
            angel_dict[__Prule1__].append([pro_seq[i:i + 3], (i+1, i + 3)])
    for combi in __Prule2__.combi:
        it = re.finditer(combi,pro_seq)
        locations = [x.start() for x in it]
        for i in locations:
            angel_dict[__Prule2__].append([pro_seq[i:i+3],(i+1,i+3)])
    angel_dict_final = {__Prule1__: [], __Prule2__: []}
    for li in angel_dict:
        for subli in angel_dict[li]:
            if subli not in angel_dict_final[li]:
                angel_dict_final[li].append(subli)
    return angel_dict_final


def D_search_rule(pro_seq):
    angel_dict = {__Drule1__:[],__Drule2__:[],__Drule3__:[]}
    for combi in __Drule1__.combi:
        it = re.finditer(combi, pro_seq)
        locations = [x.start(1) for x in it]
        for i in locations:
            angel_dict[__Drule1__].append([pro_seq[i:i+4],(i+1,i+4)])
    for combi in __Drule2__.combi:
        it = re.finditer(combi,pro_seq)
        locations = [x.start(1) for x in it]
        for i in locations:
            angel_dict[__Drule2__].append([pro_seq[i:i+4],(i+1,i+4)])
    for combi in __Drule3__.combi:
        it = re.finditer(combi,pro_seq)
        locations = [x.start(1) for x in it]
        for i in locations:
            angel_dict[__Drule3__].append([pro_seq[i:i+4],(i+1,i+4)])
    return angel_dict


def create_dat_section(seqname,sequnces_dict):
    seq = sequnces_dict[seqname]
    Dangels = D_search_rule(seq)
    Pangels = P_search_rule(seq)

    # count_angels = 1
    # Pangels_txt = ""
    # for rule in Pangels:
    #     for t in Pangels[rule]:
    #         Pangels_txt+=" " * int(5 - len(str(count_angels))) + str(count_angels)
    #         for pos in range(t[1][0], t[1][1]+1):
    #             Pangels_txt+=" " * int(5 - len(str(pos))) + str(pos)
    #         Pangels_txt+="   " + rule.angel1 + "   " + rule.strength1+"\n"
    #         count_angels += 1
    # Pangels_txt = (" "*int(12-len(str(count_angels-1)))+str(count_angels-1)+" bond angle trios."+"\n"+Pangels_txt)

    Dangels_txt = ""
    count_angels=1
    for rule in Dangels:
        for t in Dangels[rule]:
            Dangels_txt+=(" " * int(5 - len(str(count_angels))) + str(count_angels))
            for pos in range(t[1][0], t[1][1]+1):
                Dangels_txt+=(" " * int(5 - len(str(pos))) + str(pos))
            Dangels_txt+=(
                "   " + rule.angel1 + "   " + rule.strength1 + "   " + rule.angel2 + "   " + rule.strength2 + "\n")
            count_angels += 1
    Dangels_txt = " " * int(12 - len(str(count_angels - 1))) + str(count_angels - 1) + " dihedral quartets." + "\n" + Dangels_txt
    return Dangels_txt


def write_his(last_index,segment):
    last_index_number=last_index
    count = 0
    txt = ""
    for resi in segment:
        count += 1
        if resi == "H":
            last_index_number+=1
            txt += ((" " * int(8 - len(str(last_index_number)))) +
                                str(last_index_number) + (" "*int(5-len(str(count)))))
            txt+=(str(count))
            txt+=(str(3*" ")+"0.000")
            txt+=("\n")
    return txt,last_index_number


def changeAromaticToAlanin(sequnces_dict):
    for seq in sequnces_dict:
        print(sequnces_dict[seq])
        repStr = re.sub(r"[WFY]", "A", sequnces_dict[seq])
        sequnces_dict[seq]=repStr
        print(sequnces_dict[seq])
    return sequnces_dict


def create_dat(directory,eps_dirs,segment,sequnces_dict):
            dat_origin = open(directory+eps_dirs+".dat","r").read()
            hps_contacts = open(directory+"hps.dat","r").read()
            dat_write = open(directory+"new_file"+".dat","w+")
            dat_rows = dat_origin.split("\n")
            regular_sections = ["Hookean","electrostatic residues."]
            empty_sections = ["repulsive pairs."]
            unique_sections = ["contacts.","dihedral quartets.","atom positions","bond angle"]
            write_nrows = 0
            write_bonds =0
            write_electro = 0
            electro_rows = ""
            electrostatic_number=0
            his_electrostatic=""
            for row in dat_rows:
                if regular_sections[0] in row:
                    write_nrows = int(re.search("\s*([0-9]*)\s.*",row).group(1))
                    dat_write.write(row+"\n")
                    if len(regular_sections)!=1:
                        regular_sections.remove(regular_sections[0])
                elif write_nrows !=0:
                    dat_write.write(row+"\n")
                    write_nrows -=1
                elif write_bonds !=0:
                    without_end = re.match("(.*)(.{7}20\.000)",row).group(1)
                    dat_write.write(without_end+"1.745  10.000"+"\n")
                    write_bonds -= 1
                elif empty_sections[0] in row:
                    print(empty_sections[0])
                    dat_write.write(" "*11+"0 "+empty_sections[0]+"\n")
                    if len(empty_sections) != 1:
                        empty_sections.remove(empty_sections[0])
                elif "bond angle" in row:
                    write_bonds = int(re.search("\s*([0-9]*)\s.*", row).group(1))
                    dat_write.write(row+"\n")
                elif "dihedral quartets." in row:
                    dat_write.write(create_dat_section(segment,sequnces_dict))
                elif "contact" in row:
                    hps_rows = 0
                    for r in hps_contacts.split("\n"):
                        hps_rows+=1
                    hps_rows-=1
                    dat_write.write(" "*(12-len(str(hps_rows)))+str(hps_rows)+" contacts."+"\n")
                    dat_write.write(hps_contacts)
                # elif "electrostatic residues" in row:
                #     write_electro = int(re.search("\s*([0-9]*)\s.*",row).group(1))
                #     his_electrostatic,electrostatic_number = write_his(write_electro,sequnces_dict[segment])
                # elif write_electro !=0:
                #     electro_rows+=row+"\n"
                #     write_electro-=1
                elif "atom positions" in row:
                    write_nrows = int(re.search("\s*(\d*)\s.*", row).group(1))+2
                    # dat_write.write(" "*(12-len(str(electrostatic_number)))+str(electrostatic_number)+" electrostatic residues."+"\n")
                    # dat_write.write(electro_rows+his_electrostatic)
                    dat_write.write(row + "\n")
                else:
                    continue


def iterdir(directory):
    directories_toiter =[directory+"/"]
    while directories_toiter:
        dir = directories_toiter.pop()
        if glob.glob(dir+"/*.dat"):
            segment = re.search(".*(seg[ABC]).*",dir).group(1)
            epsilon = re.search(".*(eps..).*",dir).group(1)
            create_dat(dir,epsilon,segment,segment_seqs)
        else:
            for subdir in [x for x in os.walk(dir)][0][1]:
                directories_toiter.append(dir+subdir+"/")


iterdir("/home_e/dor/CPEB4/coarse_grain/HIS_charges/HIS_nocharge/")
import re

segment_seqs = {"segA":"CNTGNKSAFPVRFHPHLQPPHHHQNATPSPAAFINNNTAANGSSAGSAWLFPAPATHNIQDEILGSEKAKSQQQEQQDPLEKQQLSPSPGQEAGILPETEKAKSEENQGDNSSENGNGKEKIRIESPVLTGFDYQEATGLGTSTQPLTSSC",
                "segB":"CSSLTGFSNWSAAIAPSSSTIINEDASFFHQGGVPAASANNGALLFQNFPHHVSPGFGGSFSPQIGPLSQHHPHHPHFQHHHSQHQQQRRSPASPHPPPFTHRNAAFNQLPHLANNLNKPPSPWSSYQSPSPTPSSSWSPGGGGYGGWGGSQGRDHRRC",
                "segC":"CLNGGITPLNSISPLKKNFASNHIQLQKYARPSSAFAPKSWMEDSLNRADNIFPFPDRPRTFDMHSLESSLIDIMRAENDTIKARTYGRRRGQSSLFPMEDGFLDDGRGDQPLHSGLGSPHSFSHQNC"}

for i in segment_seqs:
    length = len(segment_seqs[i])
    aromatic = (len(re.findall("[FWY]",segment_seqs[i]))*(100/length), len(re.findall("[FWY]",segment_seqs[i])))
    proline = (len(re.findall("P",segment_seqs[i]))*(100/length),len(re.findall("P",segment_seqs[i])))
    plus = (len(re.findall("[RK]",segment_seqs[i]))*(100/length),len(re.findall("[RK]",segment_seqs[i])))
    minus = (len(re.findall("[ED]",segment_seqs[i]))*(100/length),len(re.findall("[ED]",segment_seqs[i])))
    histidine = (len(re.findall("H",segment_seqs[i]))*(100/length),len(re.findall("H",segment_seqs[i])))
    print("\n",i,"length:",length)
    print("aromatic:",aromatic)
    print("proline:",proline)
    print("plus:",plus)
    print("minus:",minus)
    print("histidine:",histidine)
    print("NCPR:",((plus[1]+histidine[1]/2-minus[1]))/length)
    print("FCR:",(aromatic[1]*0.5+plus[1]+minus[1])/length)

for i in segment_seqs:
    string = ""
    for letter in segment_seqs[i]:
        if letter == "H":
            string+="K"
        else:
            string+= letter
    print(i,string)
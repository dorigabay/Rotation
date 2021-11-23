import re

seq = "AFPPGPTTHSGPSPRPDDGPRPPPFFDVA"


class __rule1__:
    seq = ["([^P]{2}P.)", "P[^P](P[^P]P[^P])"]
    # also = "P[^P](P[^P]P[^P]"
    angel = "3.900"


class __rule2__:
    seq = ["(P[^P]P.)", "([^P]PP.)", "(PPP.)"]
    # seq2 = "([^P]PPP?)"
    # seq3 = "(PPPP?)"
    angel = "4.900"


def search_rule(pro_seq):
    rule1_search = re.findall("|".join(__rule1__.seq), pro_seq)
    rule2_search = re.findall("|".join(__rule2__.seq), pro_seq)
    angel_dict = {}
    angel_dict[__rule1__.angel] = rule1_search
    angel_dict[__rule2__.angel] = rule2_search
    print(angel_dict)
    for i in angel_dict:
        print(i)
        all_matches = []
        for t in angel_dict[i]:
            for s in t:
                print(s)
                if s=="":
                    continue
                else:
                    print(s)
                    start = re.match(s,pro_seq).start()
                    end = re.match(s,pro_seq).end()
                    all_matches.append((s,(start,end)))
        angel_dict[i] = all_matches
    return angel_dict


# def exp_dat(seq):


print(search_rule(seq))

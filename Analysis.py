import sys


if len(sys.argv) < 2:
    exit ("Not enough arguments")


arguments_given = dict(zip(sys.argv[2::2], sys.argv[3::2]))

if sys.argv[1] == "atomistic":
    directory, segment, start_step, end_step,outdir = ["" for x in range(5)]
    if "-f" not in arguments_given.keys():
        exit("""---------------Missing input directory ('-f')---------------""")
    else:
        directory = str(arguments_given["-f"]+"/")
    if "-seg" not in arguments_given.keys():
        exit("---------------Missing segment name ('-seg')---------------")
    else:
        segment = str(arguments_given["-seg"])
    if "-start" in arguments_given.keys():
        start_step = int(arguments_given["-start"])
    else: start_step = None
    if "-end" in arguments_given.keys():
        end_step = int(arguments_given["-end"])
    else: end_step = None
    if "-out" in arguments_given.keys():
        outdir = str(arguments_given["-out"])
    else: outdir = None
    from simuAnalysis2 import save_xtcAnalysis
    save_xtcAnalysis(directory, segment, outdir, start_step, end_step)

elif sys.argv[1] == "coarse":
    directory, start_step,end_step = ["" for x in range(3)]
    if "-f" not in arguments_given.keys():
        exit("""---------------Missing input directory ('-f')---------------""")
    else:
        directory = str(arguments_given["-f"]+"/")
    if "-start" in arguments_given.keys():
        start_step = int(arguments_given["-start"])
    else: start_step = None
    if "-end" in arguments_given.keys():
        end_step = int(arguments_given["-end"])
    else: end_step = None
    from simuAnalysis2 import zip2pdb, iterPDBAnalysis
    zip2pdb(directory)
    iterPDBAnalysis(directory, start_step,end_step)

elif sys.argv[1] != "coarse" or sys.argv[1] != "atomistic":
    exit("---------------Missing simulation type ('atomistic' or 'coarse')---------------")
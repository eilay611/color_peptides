import os
import sys
sys.path.insert(0, "/home/oem/PycharmProjects/lab_projects")
#import pdb_guide_Eilay_jj
spike_seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
spike_seq = spike_seq[:613]+"G"+spike_seq[614:]
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)
import platform
progect_path = os.path.dirname(os.path.abspath("coloring.py"))


def get_atom_features(pdbID):
    """
    input: pdbID && the pdb file in the same directory as this script
    without the .pdb

    output: the atom features of a pbd file
    lists of list that represents various properties of each of the atoms
    in the protein structure file. For each atom eight features are stored as follows:
    0 ATOM_serial Int The serial number of the atom in the pdb
    1 ATOM_symbol Str The name of the atom in abbreviated manner
    2 Res_number Int Residue sequence number
    3 RES_name Str Abbreviated residue name, ARG for Argenine,GLY for Glycine ..
    4 Chain_id Str The chain identifier. Some protein complexes composed of several chains
    5 X Float Orthogonal coordinates for X in Angstroms
    6 Y Float Orthogonal coordinates for Y in Angstroms
    7 Z Float Orthogonal coordinates for Z in Angstroms
    """
    features_matrix=[]#define empty tuple
    IN=open(pdbID+'.pdb','r')
    for line in IN.readlines():
        one_line_vector=()
        if line.startswith('ATOM'):
            #0123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
            atom_serial=int(line[6:11].replace(" ",""))
            atom_simbol=line[12:16].replace(" ","")

            res_number = int(line[22:26].replace(" ", ""))
            res_name=line[17:20].replace(" ","")

            chain_id = line[21:22].replace(" ", "")

            x = float(line[30:38].replace(" ", ""))
            y = float(line[38:46].replace(" ", ""))
            z = float(line[46:54].replace(" ", ""))
            one_line_vector = (atom_serial,atom_simbol,res_number,res_name,chain_id,x,y,z)  # define atom features tuple
            #print one_line_vector
            features_matrix.append(one_line_vector)
    features_matrix=tuple(features_matrix)
    IN.close()
    return features_matrix

def get_3id_to_1id_or_1id_from_3id(amino_str):
    """
    The function gets a string of amino acid name in one(‘V’) or in three letters(‘VAL’) and returns a string of an
    id of amino acids in three letters(‘VAL’) or in one letter(‘V’), respectively
    """
    dict_amino_3to1 ={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    if len(amino_str) == 3:
         return dict_amino_3to1[amino_str]
    elif len(amino_str) == 1:
        for key, value in dict_amino_3to1.items():
            if value == amino_str:
                return key
    else:
        raise ValueError("not 3 or 1 len")

def get_coordinates_of_chain(atom_features,chain_ID):
    """
    input:
     atom_features tuple and a chain ID

    output:
    dictionary of dictionaries (D_O_D) where the outer keys are the tuples (Res_number,
    amino_acid_type) and the values are dictionaries where the keys are the atoms names
    and the values are tuples with the x,y,z coordinates of the atom from the given chain.
    """
    end_dict_of_dict = {}
    for tup in atom_features:
        if tup[4] != chain_ID:
            continue
        outer_key = (tup[2],tup[3])
        inner_key = tup[1]
        if(outer_key in end_dict_of_dict):
            end_dict_of_dict[outer_key][inner_key] = (tup[5],tup[6],tup[7])
        else:
            end_dict_of_dict[outer_key]={inner_key:(tup[5],tup[6],tup[7])}
    return end_dict_of_dict

def seq(D_of_D_chain="",pdb_file_path="",chain =""):
    """
    input:
    DoD_chain1- a dictionary of dictionaries holding the atoms coordinate for each residue from chain1
    or:
    pdb_file_path - path of pdb
    chain - chain of pdb
    output:
    tupel:
    1 the seq(of orderd residued) of the DOD of the ordered_residues
    """
    if D_of_D_chain != "":
        assert chain =="" and pdb_file_path ==""
        keys = list(D_of_D_chain.keys())
        resideues = ""
        privius_place = keys[0][0]
        counter = 0
        for place , res  in keys :
            if counter==0:
                try:
                    resideues += "X"* (place-1)# fixing aligmnent on the beggining
                except TypeError:
                    resideues += "X"* (int(place[:-1])-1)# fixing 1A
                    print("i use "+ str(place[:-1])+ "insted" + str(place)  )
            jump = place - privius_place
            resideues += "X"*(jump-1) + get_3id_to_1id_or_1id_from_3id(res)
            privius_place = place
            counter+=1
    else:
        assert  chain != "" and pdb_file_path != ""
        atm = get_atom_features(pdb_file_path.rstrip(".pdb"))
        D_of_D_chain = get_coordinates_of_chain(atm,chain)
        keys = list(D_of_D_chain.keys())
        resideues = ""
        privius_place = keys[0][0]
        counter = 0
        for place , res  in keys :
            if counter==0:
                try:
                    resideues += "X"* (place-1)# fixing aligmnent on the beggining
                except TypeError:
                    resideues += "X"* (int(place[:-1])-1)# fixing aligmnent on the beggining
                    print("i use "+ str(place[0])+ "insted" + str(place)  )
            jump = place - privius_place
            resideues += "X"*(jump-1) + get_3id_to_1id_or_1id_from_3id(res)
            privius_place = place
            counter+=1

    return resideues

def run(fhandle, option_list):
    """
    by pdb_tools
    @article {Rodrigues483305,
    author = {Rodrigues, Jo{\~a}o P.G.L.M. and Teixeira, Jo{\~a}o M.C. and Trellet, Mika{\"e}l and Bonvin, Alexandre M.J.J.},
    title = {pdb-tools: a swiss army knife for molecular structures},
    elocation-id = {483305},
    year = {2018},
    doi = {10.1101/483305},
    publisher = {Cold Spring Harbor Laboratory},
    abstract = {The pdb-tools are a collection of Python scripts for working with molecular structure data in the PDB format. They allow users to edit, convert, and validate PDB files, from the command-line, in a simple but efficient manner. The pdb-tools are implemented in Python, without any external dependencies, and are freely available under the open-source Apache License at https://github.com/haddocking/pdb-tools/ and on PyPI (https://pypi.org/project/pdb-tools/).},
    URL = {https://www.biorxiv.org/content/early/2018/12/04/483305},
    eprint = {https://www.biorxiv.org/content/early/2018/12/04/483305.full.pdf},
    journal = {bioRxiv}
    }
    function from pdb_fixinsert.py

    Delete insertion codes (at specific residues).

    By default, removes ALL insertion codes on ALL residues. Also bumps
    the residue numbering of residues downstream of each insertion.

    This function is a generator.

    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.

    option_list : list
        List of insertion options to act on.
        Example ["A9", "B12"]. An empty list performs the default
        action.

    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """

    option_set = set(option_list)  # empty if option_list is empty

    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    for line in fhandle:

        if line.startswith(records):
            res_uid = line[17:27]  # resname, chain, resid, icode
            id_res = line[21] + line[22:26].strip()  # A99, B12
            has_icode = line[26].strip()  # ignore ' ' here

            # unfortunately, this is messy but not all PDB files follow a nice
            # order of ' ', 'A', 'B', ... when it comes to insertion codes..
            if prev_resi != res_uid:  # new residue
                # Does it have an insertion code
                # OR have we seen this chain + resid combination before?
                # #2 to catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
                if (has_icode or id_res in seen_ids):
                    # Do something about it
                    # if the user provided options and this residue is in them
                    # OR if the user did not provide options
                    if (not option_set) or (id_res in option_set):
                        clean_icode = True
                    else:
                        clean_icode = False
                else:
                    clean_icode = False

                prev_resi = res_uid

                if id_res in seen_ids:  # offset only if we have seen this res.
                    offset += 1

            if clean_icode:  # remove icode
                line = line[:26] + ' ' + line[27:]

            # Modify resid if necessary
            resid = int(line[22:26]) + offset
            line = line[:22] + str(resid).rjust(4) + line[26:]
            seen_ids.add(id_res)

            # Reset offset on TER
            if line.startswith('TER'):
                offset = 0

        yield line



def clean_insertion(pdb_file_path,pdb_output_path=None):
    """
    :param pdb_file_path: just the pdb file path
    :param pdb_output_path: the output file path that is pdb_file_path+".myout" if None is supply
    :return:
    """
    flag = False #is the output path is the same as the input path
    if pdb_output_path is None:
        pdb_output_path = pdb_file_path+".myout"
    elif pdb_output_path == pdb_file_path:
        pdb_output_path = pdb_file_path + ".myout"
        flag = True
    generator_of_new_pdb = run(open(pdb_file_path),[])# thanks to pdb-tools
    with open(pdb_output_path,"w") as a:
        for line in  generator_of_new_pdb:
            a.write(line)
    if flag:
        os.remove(pdb_file_path)
        os.rename(pdb_output_path,pdb_output_path.lstrip(".myout"))

def clean_insertion_folder(pdb_folder_path,pdb_output_folder_path=None):
    for file in os.listdir(pdb_folder_path):
        if file.endswith(".pdb"):
            clean_insertion(os.path.join(pdb_folder_path,file), pdb_output_path=os.path.join(pdb_output_folder_path,file))

def generate_seq_db(protein_seq1,atm_pdb_seq2,output_folder_path,muscle_folder_path,startsswith = ""):
    """
    will take tow string and make tham alighment using musscel and will return a pandas dataframe of the sequence alighment
    :param protein_seq1:
    :param atm_pdb_seq2:
    :param output_folder_path:
    :param muscle_folder_path: str the name of the folder that contain the binary file of muscle,
     the name of the binary file will be muscle, if muscle_folder_path = "" that mean that the muscle_folder_path is in
     the environmental path
    :return:
    """
    record_protein = SeqRecord(
    Seq(protein_seq1),
    id = "protein",
    name = "",
    description = "")
    record_atm = SeqRecord(
    Seq(atm_pdb_seq2),
    id = "atm_pdb",
    name = "",
    description = "")
    sequences = []
    sequences.append(record_atm)
    sequences.append(record_protein)
    with open(output_folder_path + "/input_seqs.fasta", "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")
    muscel_name = ""
    if platform.system() == "linux":
        muscel_name = "./muscle"
        sep = ";"
    elif platform.system() == "Windows":
        muscel_name = "muscle.exe"
        sep = "&"
    else:
        raise ValueError("undefinded system")
    print("cd {} {} {} -in {} -out {}".format(muscle_folder_path,sep,muscel_name,os.path.join(output_folder_path , "input_seqs.fasta"),os.path.join(output_folder_path , startsswith + "align_by_muscle_seqs.fasta")))
    os.system("cd {} {} {} -in {} -out {}".format(muscle_folder_path,sep,muscel_name,os.path.join(output_folder_path , "input_seqs.fasta"),os.path.join(output_folder_path , startsswith + "align_by_muscle_seqs.fasta")))
    seqs_db = pd.DataFrame()
    for record in SeqIO.parse(os.path.join(output_folder_path ,startsswith + "align_by_muscle_seqs.fasta"), "fasta"):
        seqs_db[record.id] = list(str(record.seq))
    seqs_db.index = seqs_db.index + 1
    seqs_db.columns = ["atm_pdb_seq","protein_seq"]
    for i in seqs_db.index:
        subseries_of_the_atm_pdb_count = seqs_db.loc[:i, "atm_pdb_seq"]  # the subseries of the atm pdb count
        pdb_atm_count = (subseries_of_the_atm_pdb_count != "-").sum()
        seqs_db.loc[i, "atm_count"] = int(pdb_atm_count)
        subseries_of_the_atm_pdb_count = seqs_db.loc[:i, "protein_seq"]  # the subseries of the atm pdb count
        pdb_atm_count = (subseries_of_the_atm_pdb_count != "-").sum()
        seqs_db.loc[i, "pep_count"] = int(pdb_atm_count)
    return seqs_db


def get_pymol_comand_for_color_list_of_residues_pic_and_movie(list_of_residues,chains,colors,name_of_selections,output_folder,together=False,take_unique_chain_together=True,len_movie=0,dpi=150,ray=0,quiet=1):
    """
    THIS IS THE SUBMAIN FUNCTION
    it color a list of residues on the pdb file path by the chains you give and the colors you give the name of selection
    will apper in the pymol
    :param list_of_residues:[[40, 41, 42, 43, 45, 46, 47, 48, 49], []]
    :param chains:['A', 'B']
    :param colors:['red', 'red']
    :param name_of_selections:['check_A', 'check_B']
    :param output_folder:
    :return:
    """
    assert len(name_of_selections) == len(list_of_residues) == len(colors) == len(name_of_selections),"the (name_of_selections,list_of_residues,colors) have to be the same length... yours is {},{},{}".format( len(name_of_selections),len(list_of_residues),len(name_of_selections))
    pymol.finish_launching()
    cmd.color("grey70","all")
    number_of_unique_chain = len(set(chains))
    print(str(number_of_unique_chain))
    counter=1
    for list_of_residue,chain,color,name_of_selection in zip(list_of_residues,chains,colors,name_of_selections):
        if list_of_residue!=[]:
            list_of_residue = [str(i) for i in list_of_residue]
            cmd.select(name_of_selection,"chain "+ chain +" and resi "+",".join(list_of_residue) )
            if type(color) is str:
                cmd.color(color, name_of_selection)
            elif type(color) is list:
                assert len(color) == 3,"if enter a list in the color the len have to be three as RGB color expected to be"
                cmd.set_color("colour_name", color)
                cmd.color("colour_name", name_of_selection)
            cmd.deselect()

        if not together==True:
            # cmd.png(os.path.join(output_folder,name_of_selection), width=0, height=0, dpi=dpi,ray=ray,quiet=quiet)
            # if len_movie!=0:
            #     cmd.mset("try", len_movie+1)
            #     util.mroll(0, len_movie, False)
            #     movie.produce(os.path.join(output_folder,name_of_selection)+".mpg")
            if counter == number_of_unique_chain or take_unique_chain_together==False:
                cmd.png(os.path.join(output_folder, "_".join(name_of_selection.split("_")[:-1])), width=0, height=0, dpi=dpi, ray=ray,
                        quiet=quiet)
                if len_movie != 0:
                    cmd.mset("try", len_movie + 1)
                    util.mroll(0, len_movie, False)
                    movie.produce(os.path.join(output_folder, name_of_selection) + ".mpg")
                counter = 0
                cmd.color("grey70", "all")
            counter += 1

    if together==True:
        cmd.png(os.path.join(output_folder,name_of_selection), width=0, height=0, dpi=dpi,ray=ray,quiet=quiet)
        if len_movie!=0:
            cmd.mset("try", len_movie+1)
            util.mroll(0, len_movie, False)
            movie.produce(os.path.join(output_folder,name_of_selection)+".mpg")
"""
cmd.extend("get_pymol_comand_for_color_list_of_residues_pic_and_movie", get_pymol_comand_for_color_list_of_residues_pic_and_movie)
get_pymol_comand_for_color_list_of_residues_pic_and_movie(
[[40, 41, 42, 43, 45, 46, 47, 48, 49], [80, 81, 82, 83, 84, 85, 86, 87, 89, 92, 91]],
['A', 'A'],
['red', 'blue'],
['check_A', 'check_B'],
"")
"""

def color_pep_on_pdb(pdb_file_path,protein_seq,chains,list_of_residue_to_color,color,name_of_selection,output_folder_path,muscle_folder_path,together=False,len_movie=0,dpi=150,ray=0,quiet=1):
    """
    THIS IS THE MAIN FUNCTION OF THIS MODULE

    it will color the list_of_residue_to_color on the pdb from your file path in pymol by the counting of the protein
    seq that you gave by doing msa by muscle.

    :param pdb_file_path: just a str file path
    :param protein_seq: the sequence of the protein that you are working on
    :param chains: the chains on the pdb that your protein is based on i.e if it is a monomer so input the monomer chain id
     ["A"] if its a dimer that build by the chains "A","B" so input ["A","B"]. the order is relevant!


    :param list_of_residue_to_color: tupel of ( list of residue to color base on the same order by your protein_seq )    or without the tupel for a alone color
    :param color:                    tupel of ( str color as pymol names of color )                                      or without the tupel for a alone color
    :param name_of_selection:        tupel of ( str name )                                                               or without the tupel for a alone color

    :param output_folder_path: a path for all the output
    :param muscle_folder_path: str the name of the folder that contain the binary file of muscle,
     the name of the binary file will be muscle, if muscle_folder_path = "" that mean that the muscle_folder_path is in
     the environmental path
    :return:
    """
    if (type(list_of_residue_to_color) is tuple) or (type(color) is tuple) or (type(name_of_selection) is tuple):
        assert (type(color) is tuple) and (type(name_of_selection) is tuple) and (type(
            list_of_residue_to_color) is tuple), "if one of those(color,name_of_selection,list_of_residue_to_color) is tupel everything should be a tupel"

    clean_insertion(pdb_file_path,os.path.join(output_folder_path,os.path.basename(pdb_file_path)))
    the_atm_seq = ""
    len_of_chain = []
    for chain in chains:
        sequence = seq(pdb_file_path =pdb_file_path,chain =chain)
        the_atm_seq +=sequence
        len_of_chain.append((chain, len(sequence)))
    seqs_db = generate_seq_db(protein_seq,the_atm_seq,output_folder_path,muscle_folder_path)#seq aligment
    seqs_db["atm_count"] = seqs_db["atm_count"].astype(float).astype(int)
    seqs_db["pep_count"] = seqs_db["pep_count"].astype(float).astype(int)
    len_up_to_now = 1
    #add chain colon
    for (chain,len_seq) in len_of_chain:
        #seqs_db.loc[(seqs_db["atm_pdb_seq"]!="-")&(seqs_db.index>=len_up_to_now)&(seqs_db.index<=len_up_to_now+len_seq),"chain"] = chain
        seqs_db.loc[(seqs_db["atm_pdb_seq"]!="-")&(seqs_db["atm_count"]>=len_up_to_now)&(seqs_db["atm_count"]<=len_up_to_now+len_seq),"chain"] = chain
        len_up_to_now += len_seq
    #reset the atm count per chain
    for chain in chains:
        seqs_db.loc[seqs_db["chain"] == chain,"atm_count"]= range(1,len(seqs_db.loc[seqs_db["chain"]==chain,"atm_count"])+1)
    #color each chain

    if type(list_of_residue_to_color) is tuple:
        dictonery_of_prams = {"list_of_residues":[],"chains":[],"colors":[],"name_of_selections":[]}
        for list_of_residue_to_color, name_of_selection , color in zip( list_of_residue_to_color, name_of_selection , color):

            dictonery_of_pram = {"list_of_residues": [], "chains": chains, "colors": [color] * len(chains),
                                 "name_of_selections": [name_of_selection + "_" + chain for chain in chains]}
            for chain in chains:
                list_of_residue_to_color_in_atm_count = list(seqs_db.loc[(seqs_db["pep_count"].isin(
                    list_of_residue_to_color)) & (seqs_db["chain"] == chain), "atm_count"])
                dictonery_of_pram["list_of_residues"].append(list_of_residue_to_color_in_atm_count)

            dictonery_of_prams["list_of_residues"] += dictonery_of_pram["list_of_residues"]
            dictonery_of_prams["name_of_selections"] += dictonery_of_pram["name_of_selections"]
            dictonery_of_prams["colors"] += dictonery_of_pram["colors"]
            dictonery_of_prams["chains"] += dictonery_of_pram["chains"]
    else:
        dictonery_of_prams = {"list_of_residues":[],"chains":chains,"colors":[color]*len(chains),"name_of_selections":[name_of_selection+"_"+chain for chain in chains]}
        for chain in chains:
            list_of_residue_to_color_in_atm_count = list(seqs_db.loc[(seqs_db["pep_count"].isin(list_of_residue_to_color)) & (seqs_db["chain"] == chain),"atm_count"])
            dictonery_of_prams["list_of_residues"].append(list_of_residue_to_color_in_atm_count)

    get_pymol_comand_for_color_list_of_residues_pic_and_movie(dictonery_of_prams["list_of_residues"],
                                                                                           dictonery_of_prams["chains"],dictonery_of_prams["colors"],
                                                                                           dictonery_of_prams["name_of_selections"], output_folder_path,
                                                              together=together,len_movie=len_movie,dpi=dpi,ray=ray,quiet=quiet)
    print("list_of_residues - "+ str(dictonery_of_prams["list_of_residues"]),"\nchains - "+str(dictonery_of_prams["chains"]),"\ncolors - "+str(dictonery_of_prams["colors"]),"\nname_of_selections - "+str(dictonery_of_prams["name_of_selections"]), "\noutput_folder_path - "+str(output_folder_path))

"""
# Make command available in PyMOL
cmd.extend("color_pep_on_pdb", color_pep_on_pdb)
#get_pymol_comand_for_color_list_of_residues_take_tow_pics_and_movie([["259","260","261"]],["A"],["red"],["checker"],"/home/oem/Downloads")
color_pep_on_pdb(
"/home/oem/Downloads/2fk0_H5N1_vietnam_H5.pdb",
"MEKIVLLFAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKANPVNDLCYPGDFNDYEELKHLLSRINHFEKIQIIPKSSWSSHEASLGVSSACPYQGESSFFRNVVWLIKKNSTYPTIKRSYNNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPRIATRSKVNGQSGRMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSTIMKSELEYGNCNTKCQTPMGAINSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNSPQRETRGLFGAIAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNECMESVRNGTYDYPQYSEEARLKREEISGVKLESIGIYQILSIYSTVASSLALAIMVAGLSLWMCSNGSLQCRICI",
["A","B"],
[40,41,42,43,45,46,47,48,49],
"red",
"check",
"/home/oem/Downloads/0",
"/home/oem/musel")
"""

"""
# Make command available in PyMOL
cmd.extend("color_pep_on_pdb", color_pep_on_pdb)
#get_pymol_comand_for_color_list_of_residues_take_tow_pics_and_movie([["259","260","261"]],["A"],["red"],["checker"],"/home/oem/Downloads")
color_pep_on_pdb(
"/home/oem/Downloads/2fk0_H5N1_vietnam_H5.pdb",
"MEKIVLLFAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKANPVNDLCYPGDFNDYEELKHLLSRINHFEKIQIIPKSSWSSHEASLGVSSACPYQGESSFFRNVVWLIKKNSTYPTIKRSYNNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPRIATRSKVNGQSGRMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSTIMKSELEYGNCNTKCQTPMGAINSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNSPQRETRGLFGAIAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNECMESVRNGTYDYPQYSEEARLKREEISGVKLESIGIYQILSIYSTVASSLALAIMVAGLSLWMCSNGSLQCRICI",
["A","B"],
([40,41,42,43,45,46,47,48,49],[80, 81, 82, 83, 84, 85, 86, 87, 89, 92, 91]),
("red","blue"),
("check","check2"),
"/home/oem/Downloads/0",
"/home/oem/musel")
"""
dict_of_comands={'H1': 'color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/3lzg_H1N1_california2009_H1.pdb","MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI",[\'A\', \'B\'],([81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565], [41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565]),(\'red\', \'blue\', \'green\', \'magenta\', \'yellow\', \'cyan\'),(\'H1_live\', \'H1_dead\', \'H1_pbslive\', \'H1_rapalive\', \'H1_rapadead\'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel") ', 'H3': 'color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/3vun_H3N2_x31_H3.pdb","MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI",[\'A\', \'B\'],([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506], [157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481], [336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500], [221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440]),(\'red\', \'blue\', \'green\', \'magenta\', \'yellow\', \'cyan\'),(\'H3_live\', \'H3_dead\', \'H3_pbslive\', \'H3_rapalive\', \'H3_rapadead\'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel") ', 'H5': 'color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/2fk0_H5N1_vietnam_H5.pdb","MEKIVLLFAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKANPVNDLCYPGDFNDYEELKHLLSRINHFEKIQIIPKSSWSSHEASLGVSSACPYQGESSFFRNVVWLIKKNSTYPTIKRSYNNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPRIATRSKVNGQSGRMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSTIMKSELEYGNCNTKCQTPMGAINSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNSPQRETRGLFGAIAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNECMESVRNGTYDYPQYSEEARLKREEISGVKLESIGIYQILSIYSTVASSLALAIMVAGLSLWMCSNGSLQCRICI",[\'A\', \'B\'],([156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565], [146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195], [306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325], [531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500]),(\'red\', \'blue\', \'green\', \'magenta\', \'yellow\', \'cyan\'),(\'H5_live\', \'H5_dead\', \'H5_pbslive\', \'H5_pbsdead\', \'H5_rapalive\', \'H5_rapadead\'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel") ', 'Cal_N1': 'color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/3nss_H1N1_california2009_N1.pdb","MNPNQKIITIGSVCMTIGMANLILQIGNIISIWISHSIQLGNQNQIETCNQSVITYENNTWVNQTYVNISNTNFAAGQSVVSVKLAGNSSLCPVSGWAIYSKDNSVRIGSKGDVFVIREPFISCSPLECRTFFLTQGALLNDKHSNGTIKDRSPYRTLMSCPIGEVPSPYNSRFESVAWSASACHDGINWLTIGISGPDNGAVAVLKYNGIITDTIKSWRNNILRTQESECACVNGSCFTVMTDGPSNGQASYKIFRIEKGKIVKSVEMNAPNYHYEECSCYPDSSEITCVCRDNWHGSNRPWVSFNQNLEYQIGYICSGIFGDNPRPNDKTGSCGPVSSNGANGVKGFSFKYGNGVWIGRTKSISSRNGFEMIWDPNGWTGTDNNFSIKQDIVGINEWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKENTIWTSGSSISFCGVNSDTVGWSWPDGAELPFTIDK",[\'A\'],([21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185], [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440]),(\'red\', \'blue\', \'green\', \'magenta\', \'yellow\', \'cyan\'),(\'Cal_N1_live\', \'Cal_N1_pbslive\', \'Cal_N1_rapalive\', \'Cal_N1_rapadead\'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel") ', 'N2': 'color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/","MNPNQKIITIGSVSLTIATVCFLMQIAILVTTVTLHFKQYECDSPASNQVMPCEPIIIERNITEIVYLNNTTIEKEICPKVVEYRNWSKPQCQITGFAPFSKDNSIRLSAGGDIWVTREPYVSCDHGKCYQFALGQGTTLDNKHSNDTIHDRIPHRTLLMNELGVPFHLGTRQVCIAWSSSSCHDGKAWLHVCITGDDKNATASFIYDGRLVDSIGSWSQNILRTQESECVCINGTCTVVMTDGSASGRADTRILFIEEGKIVHISPLSGSAQHVEECSCYPRYPGVRCICRDNWKGSNRPVVDINMEDYSIDSSYVCSGLVGDTPRNDDRSSNSNCRNPNNERGNQGVKGWAFDNGDDVWMGRTISKDLRSGYETFKVIGGWSTPNSKSQINRQVIVDSDNRSGYSGIFSVEGKSCINRCFYVELIRGRKQETRVWWTSNSIVVFCGTSGTYGTGSWPDGANINFMPI",[\'A\'],([82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461], [312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331], [96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290], [311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330], [81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400], [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445]),(\'red\', \'blue\', \'green\', \'magenta\', \'yellow\', \'cyan\'),(\'N2_live\', \'N2_dead\', \'N2_pbslive\', \'N2_pbsdead\', \'N2_rapalive\', \'N2_rapadead\'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel") ', 'Vie_N1': 'color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/2hty_H5N1_vietnam_N1.pdb","MNPNQKIITIGSICMVTGIVSLMLQIGNMISIWVSHSIHTGNQHQSEPISNTNFLTEKAVASVKLAGNSSLCPINGWAVYSKDNSIRIGSKGDVFVIREPFISCSHLECRTFFLTQGALLNDKHSNGTVKDRSPHRTLMSCPVGEAPSPYNSRFESVAWSASACHDGTSWLTIGISGPDNGAVAVLKYNGIITDTIKSWRNNILRTQESECACVNGSCFTVMTDGPSNGQASHKIFKMEKGKVVKSVELDAPNYHYEECSCYPNAGEITCVCRDNWHGSNRPWVSFNQNLEYQIGYICSGVFGDNPRPNDGTGSCGPVSSNGAYGVKGFSFKYGNGVWIGRTKSTNSRSGFEMIWDPNGWTETDSSFSVKQDIVAITDWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKESTIWTSGSSISFCGVNSDTVGWSWPDGAELPFTIDK",[\'A\'],([16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380], [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355], [61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430]),(\'red\', \'blue\', \'green\', \'magenta\', \'yellow\', \'cyan\'),(\'Vie_N1_live\', \'Vie_N1_dead\', \'Vie_N1_pbslive\', \'Vie_N1_rapalive\', \'Vie_N1_rapadead\'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel") '}
cmd.extend("color_pep_on_pdb", color_pep_on_pdb)
# cmd.set_view ((
#      0.136094809,    0.584857523,   -0.799637854,
#     -0.843417943,   -0.355043322,   -0.403225690,
#     -0.519735098,    0.729304910,    0.444958538,
#     -0.000006974,    0.000034586, -409.834259033,
#     12.266601562,   41.237243652,  -20.679107666,
#    330.361083984,  489.309509277,  -20.000000000 ))
#H1

#color_pep_on_pdb("C:\\Users\\Public\\Downloads\\3lzg.pdb","MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI",['A', 'B'],([81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565], [41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565]),('red', 'blue', 'green', 'magenta', 'yellow', 'cyan'),('H1_live', 'H1_dead', 'H1_pbslive', 'H1_rapalive', 'H1_rapadead'),"C:\\Users\\Public\\Downloads\\check","C:\\Users\\Public\\Downloads")
color_pep_on_pdb(pdb_file_path=os.path.join(progect_path,"example","3lzg.pdb"),
                 protein_seq="MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETPSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNIPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI",
                 chains=['A', 'B'],
                 list_of_residue_to_color=([36,37,38,39,41,40,43,44],[100,101,102,103,104,105,106],[150,151,152,153,154,155],[180,181,182,183,184,185],[200,201,202,203,204,205]),
                 color=('red', 'blue', 'green', 'magenta', 'yellow'),
                 name_of_selection=('H1_red', 'H1_blue', 'H1_green', 'H1_magenta', 'H1_yellow'),
                 output_folder_path=os.path.join(progect_path,"example","outputs"),
                 muscle_folder_path=os.path.join(progect_path,"example"),
                 together=False,
                 len_movie=0,
                 dpi=150,
                 ray=0,
                 quiet=1)
# #H3
#color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/3vun_H3N2_x31_H3.pdb","MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI",['A', 'B'],([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506], [157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481], [336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500], [221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440]),('red', 'blue', 'green', 'magenta', 'yellow', 'cyan'),('H3_live', 'H3_dead', 'H3_pbslive', 'H3_rapalive', 'H3_rapadead'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel")
# #H5
#color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/2fk0_H5N1_vietnam_H5.pdb","MEKIVLLFAIVSLVKSDQICIGYHANNSTEQVDTIMEKNVTVTHAQDILEKKHNGKLCDLDGVKPLILRDCSVAGWLLGNPMCDEFINVPEWSYIVEKANPVNDLCYPGDFNDYEELKHLLSRINHFEKIQIIPKSSWSSHEASLGVSSACPYQGESSFFRNVVWLIKKNSTYPTIKRSYNNTNQEDLLVLWGIHHPNDAAEQTKLYQNPTTYISVGTSTLNQRLVPRIATRSKVNGQSGRMEFFWTILKPNDAINFESNGNFIAPEYAYKIVKKGDSTIMKSELEYGNCNTKCQTPMGAINSSMPFHNIHPLTIGECPKYVKSNRLVLATGLRNSPQRETRGLFGAIAGFIEGGWQGMVDGWYGYHHSNEQGSGYAADKESTQKAIDGVTNKVNSIIDKMNTQFEAVGREFNNLERRIENLNKKMEDGFLDVWTYNAELLVLMENERTLDFHDSNVKNLYDKVRLQLRDNAKELGNGCFEFYHKCDNECMESVRNGTYDYPQYSEEARLKREEISGVKLESIGIYQILSIYSTVASSLALAIMVAGLSLWMCSNGSLQCRICI",['A', 'B'],([156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565], [146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195], [306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325], [531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550], [51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500]),('red', 'blue', 'green', 'magenta', 'yellow', 'cyan'),('H5_live', 'H5_dead', 'H5_pbslive', 'H5_pbsdead', 'H5_rapalive', 'H5_rapadead'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel")
# #Cal_N1
#color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/3nss_H1N1_california2009_N1.pdb","MNPNQKIITIGSVCMTIGMANLILQIGNIISIWISHSIQLGNQNQIETCNQSVITYENNTWVNQTYVNISNTNFAAGQSVVSVKLAGNSSLCPVSGWAIYSKDNSVRIGSKGDVFVIREPFISCSPLECRTFFLTQGALLNDKHSNGTIKDRSPYRTLMSCPIGEVPSPYNSRFESVAWSASACHDGINWLTIGISGPDNGAVAVLKYNGIITDTIKSWRNNILRTQESECACVNGSCFTVMTDGPSNGQASYKIFRIEKGKIVKSVEMNAPNYHYEECSCYPDSSEITCVCRDNWHGSNRPWVSFNQNLEYQIGYICSGIFGDNPRPNDKTGSCGPVSSNGANGVKGFSFKYGNGVWIGRTKSISSRNGFEMIWDPNGWTGTDNNFSIKQDIVGINEWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKENTIWTSGSSISFCGVNSDTVGWSWPDGAELPFTIDK",['A'],([21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185], [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440]),('red', 'blue', 'green', 'magenta', 'yellow', 'cyan'),('Cal_N1_live', 'Cal_N1_pbslive', 'Cal_N1_rapalive', 'Cal_N1_rapadead'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel")
# #N2
#color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/2hty_H5N1_vietnam_N1.pdb","MNPNQKIITIGSVSLTIATVCFLMQIAILVTTVTLHFKQYECDSPASNQVMPCEPIIIERNITEIVYLNNTTIEKEICPKVVEYRNWSKPQCQITGFAPFSKDNSIRLSAGGDIWVTREPYVSCDHGKCYQFALGQGTTLDNKHSNDTIHDRIPHRTLLMNELGVPFHLGTRQVCIAWSSSSCHDGKAWLHVCITGDDKNATASFIYDGRLVDSIGSWSQNILRTQESECVCINGTCTVVMTDGSASGRADTRILFIEEGKIVHISPLSGSAQHVEECSCYPRYPGVRCICRDNWKGSNRPVVDINMEDYSIDSSYVCSGLVGDTPRNDDRSSNSNCRNPNNERGNQGVKGWAFDNGDDVWMGRTISKDLRSGYETFKVIGGWSTPNSKSQINRQVIVDSDNRSGYSGIFSVEGKSCINRCFYVELIRGRKQETRVWWTSNSIVVFCGTSGTYGTGSWPDGANINFMPI",['A'],([82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461], [312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331], [96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290], [311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330], [81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400], [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445]),('red', 'blue', 'green', 'magenta', 'yellow', 'cyan'),('N2_live', 'N2_dead', 'N2_pbslive', 'N2_pbsdead', 'N2_rapalive', 'N2_rapadead'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel")
# #Vie_N1
#color_pep_on_pdb("/home/oem/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions/2hty_H5N1_vietnam_N1.pdb","MNPNQKIITIGSICMVTGIVSLMLQIGNMISIWVSHSIHTGNQHQSEPISNTNFLTEKAVASVKLAGNSSLCPINGWAVYSKDNSIRIGSKGDVFVIREPFISCSHLECRTFFLTQGALLNDKHSNGTVKDRSPHRTLMSCPVGEAPSPYNSRFESVAWSASACHDGTSWLTIGISGPDNGAVAVLKYNGIITDTIKSWRNNILRTQESECACVNGSCFTVMTDGPSNGQASHKIFKMEKGKVVKSVELDAPNYHYEECSCYPNAGEITCVCRDNWHGSNRPWVSFNQNLEYQIGYICSGVFGDNPRPNDGTGSCGPVSSNGAYGVKGFSFKYGNGVWIGRTKSTNSRSGFEMIWDPNGWTETDSSFSVKQDIVAITDWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKESTIWTSGSSISFCGVNSDTVGWSWPDGAELPFTIDK",['A'],([16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380], [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355], [61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430]),('red', 'blue', 'green', 'magenta', 'yellow', 'cyan'),('Vie_N1_live', 'Vie_N1_dead', 'Vie_N1_pbslive', 'Vie_N1_rapalive', 'Vie_N1_rapadead'),"/home/oem/PycharmProjects/lab_projects/rapa/seq/significant/figs","/home/oem/musel")
#





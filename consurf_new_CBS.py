#run in pymol to generate a bfactor scale from 1-9 colorized in color blinds colors

# Define a Python subroutine to colour atoms by B-factor, using predefined intervals
#import matplotlib
import os

from matplotlib import cm
import pandas as pd


def get_atom_features_db(pdb_file_path):
    """
    input:  the pdb file path
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
    IN=open(pdb_file_path,'r')
    for line in IN.readlines():
        one_line_vector=()
        if line.startswith('ATOM'):
            #0123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
            atom_serial=int(line[6:11].replace(" ",""))
            atom_simbol=line[12:16].replace(" ","")

            res_name=line[17:20].replace(" ","")
            chain_id = line[21:22].replace(" ", "")

            res_number = int(line[22:26].replace(" ", ""))

            x = float(line[30:38].replace(" ", ""))
            y = float(line[38:46].replace(" ", ""))
            z = float(line[46:54].replace(" ", ""))

            condormation_rate = float(line[55:61].replace(" ", ""))
            try:
                b_factor =  "{:.2f}".format(float(line[61:67].replace(" ", "")))
            except ValueError:
                b_factor = None
            one_line_vector = (atom_serial,atom_simbol,res_name,chain_id,res_number,x,y,z,condormation_rate,b_factor)  # define atom features tuple
            #print one_line_vector
            features_matrix.append(one_line_vector)
    features_matrix=tuple(features_matrix)
    IN.close()
    return pd.DataFrame(features_matrix,columns = ["ATOM_serial","ATOM_symbol","RES_name","Res_number","Chain_id","X","Y","Z","condormation_rate","b_factor"])

def colour_consurf_CBS(pdb_file_path=None,selection="all",minimum = None,maximum = None,n_colours = 26,colours = "Reds",output_folder_path = ""):
    """
    :param pdb_file_path:
    :param selection:
    :param minimum:
    :param maximum:
    :param n_colours:
    :param colours:
    :return:
    """

    """
    for color blinds 
    colours = [[0.058823529, 0.352941176, 0.137254902],[0.352941176, 0.68627451, 0.37254902],[0.647058824, 0.862745098, 0.62745098],[0.843137255, 0.941176471, 0.823529412],[1, 1, 1],[0.901960784, 0.823529412, 0.901960784],[0.764705882, 0.647058824, 0.803921569],[0.607843137, 0.431372549, 0.666666667],[0.470588235, 0.156862745, 0.509803922]]
    for not
        colours = [ [0.039215686, 0.490196078, 0.509803922],[0.294117647, 0.68627451, 0.745098039],[0.647058824, 0.862745098, 0.901960784],[0.843137255, 0.941176471, 0.941176471],[1, 1, 1],[0.980392157, 0.921568627, 0.960784314],[0.980392157, 0.784313725, 0.862745098],[0.941176471, 0.490196078, 0.666666667],[0.62745098, 0.156862745, 0.37254902]]
    """
    pymol.finish_launching()

    cmd.bg_color("white")
    if colours =="CBS":#color blinds color
        colours = [[0.058823529, 0.352941176, 0.137254902],[0.352941176, 0.68627451, 0.37254902],[0.647058824, 0.862745098, 0.62745098],[0.843137255, 0.941176471, 0.823529412],[1, 1, 1],[0.901960784, 0.823529412, 0.901960784],[0.764705882, 0.647058824, 0.803921569],[0.607843137, 0.431372549, 0.666666667],[0.470588235, 0.156862745, 0.509803922]]
    elif colours == "NCBS":#not color blinds color
        colours = [ [0.039215686, 0.490196078, 0.509803922],[0.294117647, 0.68627451, 0.745098039],[0.647058824, 0.862745098, 0.901960784],[0.843137255, 0.941176471, 0.941176471],[1, 1, 1],[0.980392157, 0.921568627, 0.960784314],[0.980392157, 0.784313725, 0.862745098],[0.941176471, 0.490196078, 0.666666667],[0.62745098, 0.156862745, 0.37254902]]
    elif type(colours) is list:#Colours are calculated by dividing the RGB colours by 255
        colours = colours
    elif type(colours) is str: #cmap name
        cmap = cm.get_cmap(colours, n_colours)
        colours = []
        for i in range(cmap.N):
            rgb = cmap(i)[:3]
            colours.append(rgb)
    # Colour other chains gray, while maintaining
    cmd.color("gray", selection)
    #cmd.util.cnc()    # oxygen in red, nitrogen in blue and hydrogen in white

    if pdb_file_path is None:
        # These are constants
        minimum = 0
        maximum = len(colours)
    else:
        atm = get_atom_features_db(pdb_file_path)
        minimum = atm["b_factor"].astype(float).min() #minimum of b-color
        maximum = atm["b_factor"].astype(float).max() #maximum of b-color

    n_colours = len(colours)
    # Colours are calculated by dividing the RGB colours by 255
    # RGB = [[[27,120,55],[90,174,97],[166,219,160],[217,240,211],[255,255,255],
    #        [231,212,232],[194,165,207],[153,112,171],[118,42,131]]

    bin_size = (maximum - minimum) / n_colours
    
    # Loop through colour intervals
    for i in range(n_colours):
        lower = minimum + (i + 1) * bin_size
        upper = lower + bin_size
        colour = colours[i]
        
        # Print out B-factor limits and the colour for this group
        print(lower, " - ", upper, " = ", colour)
        
        # Define a unique name for the atoms which fall into this group
        group = selection + "_group_" + str(i + 1)
        
        # Compose a selection command which will select all atoms which are
        #	a) in the original selection, AND
        #	b) have B factor in range lower <= b < upper
        sel_string = selection + " & ! b < " + str(lower)
        
        if(i < n_colours):
            sel_string += " & b < " + str(upper)
        else:
            sel_string += " & ! b > " + str(upper)
        
        # Select the atoms
        cmd.select(group, sel_string)
        
        # Create a new colour
        colour_name = "colour_" + str(i + 1)
        cmd.set_color(colour_name, colour)
        
        # Colour them
        cmd.color(colour_name, group)
    
    
    # Create new colour for insufficient sequences
    # RGB_colour = [255,255,150]
    insuf_colour = [1, 1, 0.588235294]
    cmd.set_color("insufficient_colour", insuf_colour)
    
    # Colour atoms with B-factor of 10 using the new colour
    cmd.select("insufficient", selection + " & b = 10")
    cmd.color("insufficient_colour", "insufficient")


    cmd.show("surface")
    cmd.refresh()

    if output_folder_path != "":
        cmd.set_view((
            0.959461868, 0.019938128, -0.281131893,
            0.281554550, -0.112636462, 0.952910602,
            -0.012666089, -0.993435919, -0.113684900,
            0.000380397, 0.000116260, -524.335754395,
            224.599838257, 220.184539795, 222.841293335,
            -39766.035156250, 40814.699218750, -20.000000000))

        name_of_selection = os.path.basename(pdb_file_path).replace(".pdb","")
        cmd.png(os.path.join(output_folder_path, name_of_selection+".png"), width=0, height=0, dpi=300, ray=1, quiet=1)
        cmd.mset("try", 121)
        # util.mroll(0, 120, False)
        # movie.produce(os.path.join(output_folder_path, name_of_selection + ".mpg"))


# Make command available in PyMOL
cmd.extend("colour_consurf_CBS", colour_consurf_CBS)
#run /home/oem/PycharmProjects/color_peptide_on_pdb/consurf_new_CBS.py
#colour_consurf_CBS(pdb_file_path ="/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new/unregardent/analyzide_data/1_6VYB.pdb",colours ="Reds",output_folder_path = "/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new/unregardent/analyzide_data")
colour_consurf_CBS()
# Make all groups unselected
cmd.deselect()

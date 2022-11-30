from pymol import cmd
import os
from pymol import stored

#cmd.set("surface_quality",2)
def represent_bfactorize(the_antigen=""):
    mapper={"cluster1": "6XEYKJ",
    "cluster2": "7JW0FG",
    "cluster3": "7KMHBA",
    "cluster4": "7LXYON",
    "cluster5": "7WOAGF",
    "cluster6": "7R40FD",
    "cluster7": "7PQZBA",
    "cluster8": "7BNVLH",
    "cluster9": "7DZXLH",
    "cluster10": "7AKDLH",
    "cluster11": "7RQ6KJ",
    "cluster12": "7FAELH",
    "cluster13": "7U2DLH",
    "cluster14": "8D6ZDC",
    "cluster15": "7RBVLH",
    "cluster16": "7S0ELH",
    "cluster17": "7RAQLH",
    "cluster18": "7SOBBC",
    "cluster19": "7UARGF",
    "cluster20": "7WK9ba"}
    mapper = {v:k for k,v in mapper.items()}
    cmd.set_color( "cluster1", [0.8941176470588236, 0.10196078431372549, 0.10980392156862745])
    cmd.set_color(  "cluster2", [0.21568627450980393, 0.49411764705882355, 0.7215686274509804])
    cmd.set_color(  "cluster3", [0.30196078431372547, 0.6862745098039216, 0.2901960784313726])
    cmd.set_color(  "cluster4", [0.596078431372549, 0.3058823529411765, 0.6392156862745098])
    cmd.set_color(  "cluster5", [1.0, 0.4980392156862745, 0.0])
    cmd.set_color(  "cluster6", [1.0, 1.0, 0.2])
    cmd.set_color(  "cluster7", [0.6509803921568628, 0.33725490196078434, 0.1568627450980392])
    cmd.set_color(  "cluster8", [0.9686274509803922, 0.5058823529411764, 0.7490196078431373])
    cmd.set_color(  "cluster9", [0.10196078431372549, 0.8196078431372549, 1.0])
    cmd.set_color(  "cluster10", [0.3, 0.3, 0.3])
    cmd.set_color(  "cluster11",[0.4, 0.4, 0.4])
    cmd.set_color(  "cluster12", [0.8, 0.6, 1.0])
    cmd.set_color(  "cluster13",[0.7058823529411765, 0.6274509803921569, 0.39215686274509803])
    cmd.set_color(  "cluster14", [0.7019607843137254, 0.0, 0.5254901960784314])
    cmd.set_color(  "cluster15",[0.0, 0.2, 0.8])
    cmd.set_color(  "cluster16", [0.5,0.5,0.5])
    cmd.set_color(  "cluster17",[0.6,0.6,0.6])
    cmd.set_color(  "cluster18", [0.7,0.7,0.7])
    cmd.set_color(  "cluster19",[0.8,0.8,0.8])
    cmd.set_color(  "cluster20",[0.9,0.9,0.9])
    cmd.set_color("bordo", [165, 0, 0])

    for obj in cmd.get_object_list('all'):
        try:
            mapper[obj.split("_")[0]]
        except:
            cmd.delete(obj)
            continue
        antigen = "obj "+obj+" and chain " +obj.split("_")[-1][0]
        cmd.color("blue",antigen)
        cmd.show_as("surface",antigen)
        cmd.spectrum("b", "gray75_yelloworange_bordo",antigen,None,None,1)

        if the_antigen=="NTD":
            wheats = "resi 14-33+55-83+85-87+104-115+132-139+157-167+232-239+249-277+288-304"
        elif the_antigen=="RBD":
            wheats = "resi 333-360+436-500"
        if the_antigen != "":
            cmd.color("wheat",antigen +" and {} and b<0.0001".format(wheats))
            cmd.color("palecyan",antigen +" and (not {}) and b<0.0001".format(wheats))

        cmd.color(mapper[obj.split("_")[0]],obj+" and chain "+obj.split("_")[0][4]+"+"+obj.split("_")[0][5])

        cmd.set_name(obj,mapper[obj.split("_")[0]]+"_"+obj.split("_")[0]+"_"+obj.split("_")[-1])
cmd.extend("represent_bfactorize", represent_bfactorize)

def remove_over_110_abs():
    for obj in cmd.get_object_list('all'):
        cmd.remove("obj "+obj+" and resi 110-600 and chain "+obj.split("_")[1][-2]+"+"+obj.split("_")[1][-1])
cmd.extend("remove_over_110_abs", remove_over_110_abs)
#remove_over_110_abs()
#represent_bfactorize(the_antigen="NTD")
def extract_spike():
    for obj in cmd.get_object_list('all'):
        cmd.extract("{}_ab".format(obj.split("_")[0]),"obj {} and chain {}".format(obj,obj.split("_")[1][-1]+"+"+obj.split("_")[1][-2]))
        cmd.group(obj.split("_")[0]," ".join([obj,"{}_ab".format(obj.split("_")[0])]))
cmd.extend("extract_spike", extract_spike)


def binerize(the_antigen="",cutoff=179/2,color =(210, 136, 162)):
    cmd.set_color("matte_pink", color)
    if the_antigen=="NTD":
        wheats = "resi 14-33+55-83+85-87+104-115+132-139+157-167+232-239+249-277+288-304"
    elif the_antigen=="RBD":
        wheats = "resi 333-360+436-500"
    for obj in cmd.get_object_list('all'):
        antigen = "obj "+obj+" and chain " +obj.split("_")[-1][0]
        if the_antigen != "":
            cmd.color("wheat", antigen + " and {} and b<{}".format(wheats,cutoff))
            cmd.color("palecyan", antigen + " and (not {}) and b<{}".format(wheats,cutoff))
    cmd.color("matte_pink","b>{}".format(cutoff))
cmd.extend("binerize", binerize)
#binerize(the_antigen="RBD",cutoff=int(179/2),color =(210, 136, 162))






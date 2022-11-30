#for a in cmd.get_object_list('all'):
import os
from pymol import cmd
import pickle as pkl

#####
def picture_alot_of_all(output_folder_path,name=""):
    cmd.bg_color("white")
    cmd.png(os.path.join(output_folder_path,"1_{}.png".format(name)), width=0, height=0, dpi=200, ray=0, quiet=1)
    cmd.turn("y",120)
    cmd.png(os.path.join(output_folder_path,"2_{}.png".format(name)), width=0, height=0, dpi=200, ray=0, quiet=1)
    cmd.turn("y",120)
    cmd.png(os.path.join(output_folder_path,"3_{}.png".format(name)), width=0, height=0, dpi=200, ray=0, quiet=1)
    cmd.turn("y",120)
    #back to first place
    cmd.turn("x",90)
    cmd.png(os.path.join(output_folder_path,"4_{}.png".format(name)), width=0, height=0, dpi=200, ray=0, quiet=1)

    cmd.turn("x",-45)
    cmd.png(os.path.join(output_folder_path,"5_{}.png".format(name)), width=0, height=0, dpi=200, ray=0, quiet=1)
    cmd.turn("y",120)
    cmd.png("/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new/unregardent/analyzide_data/new_model/6_{}.png".format(name), width=0, height=0, dpi=200, ray=0, quiet=1)
    cmd.turn("y",120)
    cmd.png("/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new/unregardent/analyzide_data/new_model/7_{}.png".format(name), width=0, height=0, dpi=200, ray=0, quiet=1)

    cmd.turn("y",120)
    cmd.turn("x",-45)
    #back to first place

cmd.extend("picture", picture_alot_of_all)
op="/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/cov2_280822/regardent/analyzide_data/midi_bigger_cmpnd/pdb_abs_align/representive_1percluster_pics/RBD"
op="/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/cov2_280822/regardent/analyzide_data/midi_bigger_cmpnd/pdb_abs_align/representive_1percluster_pics/RBD/b_factorize/binary"
op="/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/cov2_280822/regardent/analyzide_data/midi_bigger_cmpnd/pdb_abs_align/representive_togther_pics/NTD"
def picture_roll(output_folder_path=op,bg="white",axis = "y",number_of_shots=8,angle=None,name="",width=0, height=0, dpi=300, ray=0, quiet=1):
    cmd.bg_color(bg)
    counter = 0
    for num_pic in range(number_of_shots):
        if angle == None:
            angle=360/number_of_shots
        cmd.png(os.path.join(output_folder_path, "{}_{}.png".format(counter,name)), width=width, height=height, dpi=dpi, ray=ray, quiet=quiet)
        cmd.turn(axis, angle)
        counter+=angle

cmd.extend("picture_roll", picture_roll)

# objects = []
# with (open("/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new/unregardent/analyzide_data/strong_matrix_dict_for_pymol.pkl", "rb")) as openfile:
#     while True:
#         try:
#             objects.append(pkl.load(openfile))
#         except EOFError:
#             break
# rgb_dict = objects[0]
def color_every_resi_by_list_of_rgb_colors(rgb_dict, sub_group="all",color_to_ignore=(0.7,0.7,0.7)):
    """
    :param rgb_dict:{<chain>_<resi>:(R,G,B),...}
    :param sub_group: inner selection in your pymol ssesiion
    :param color_to_ignore: ignore this color and dont change it
    :return:
    """
    print("start")
    for key,color in rgb_dict.items():
        if color != color_to_ignore:
            chain,resi = key.split("_")
            cmd.set_color( "temp"+key, color )
            print("cmd.set_color("+ "temp"+str(chain)+"," +str(color) +")")
            print("cmd.color("+"temp"+str(chain)+","+ str(sub_group) +" and resi "+str(resi)+" and chain "+str(chain)+",1)")
            cmd.color("temp"+key, sub_group +" and resi "+resi+" and chain "+chain,1)
cmd.extend("color_every_resi_by_list_of_rgb_colors", picture_roll)

# run /home/oem/PycharmProjects/color_peptides/pic.py
# picture_roll(output_folder_path=op,axis = "y",number_of_shots=5,angle=120/5,name="",width=0, height=0, dpi=300, ray=1, quiet=1)
#picture_roll(output_folder_path=os.path.join(op,"5"),name="cluster5")
def movie_eilay(len_movie=360,output_folder_path=op,bg="white",name=""):
    cmd.bg_color(bg)
    cmd.mset("try", len_movie + 1)
    util.mroll(0, len_movie, False)
    movie.produce(os.path.join(output_folder_path, name) + ".mpg")
cmd.extend("movie_eilay", movie_eilay)


def picture_roll_and_movie(output_folder_path=op,bg="white",axis = "y",len_movie=360,number_of_shots=8,name="",width=0, height=0, dpi=300, ray=0, quiet=1):
    picture_roll(output_folder_path=output_folder_path, bg=bg, axis=axis, number_of_shots=number_of_shots, angle=None, name=name, width=width, height=height,
                 dpi=dpi, ray=ray, quiet=quiet)
    movie_eilay(len_movie=len_movie, output_folder_path=output_folder_path, bg=bg, name=name)
cmd.extend("picture_roll_and_movie", picture_roll_and_movie)


def pic_every_obgect_alone(output_folder_path=op,bg="white",name="",width=0, height=0, dpi=300, ray=0, quiet=1,shots=[],object_map_name={}):
    """
    eilay if you have groups enable tham and # the #cmd.disable("all") only the groups not the obj inside tham
    :param output_folder_path:
    :param bg:
    :param name:
    :param width:
    :param height:
    :param dpi:
    :param ray:
    :param quiet:
    :param shots:
    :param object_map_name:
    :return:
    """
    if type(shots)== []:
        shots.append(cmd.get_view())
    cmd.bg_color(bg)
    all_obgects = cmd.get_names('objects', 0, '(all)')
    cmd.disable(" ".join(all_obgects))
    counter=0
    if type(shots) == type(8):
        shots_new=[]
        for i in range (shots):
            shots_new.append(cmd.get_view())
            cmd.turn("y", 360/shots)
        shots = shots_new
    for shot in shots:
        cmd.set_view(shot)
        for obj in all_obgects:
            map_obj=obj
            if object_map_name != {}:
                try:
                    map_obj = object_map_name[obj[:6]]########the :6 is cutomize
                except:
                    continue
            cmd.enable(obj)
            cmd.png(os.path.join(output_folder_path, "{}_{}_{}.png".format(counter,name,map_obj)), width=width, height=height, dpi=dpi, ray=ray, quiet=quiet)
            cmd.disable(obj)
        counter+=1
cmd.extend("pic_every_obgect_alone", pic_every_obgect_alone)

def movie_every_obgect_alone(output_folder_path=op,bg="white",name="",len_movie=240,object_map_name={}):
    cmd.bg_color(bg)
    all_obgects = cmd.get_names('objects', 0, '(all)')
    cmd.disable(" ".join(all_obgects))
    for obj in all_obgects:
        map_obj = obj
        if object_map_name != {}:
            try:
                map_obj = object_map_name[obj[:6]]  ########the :6 is cutomize
            except:
                continue
        cmd.enable(obj)
        movie_eilay(len_movie=len_movie, output_folder_path=op, bg=bg, name=map_obj+"_"+name)
        cmd.disable(obj)
cmd.extend("movie_every_obgect_alone", movie_every_obgect_alone)

def get_groups():
    groups_names = set(cmd.get_names("objects")).difference(cmd.get_object_list("all"))
    dict_of_groups = {}
    for groups_name in groups_names:
        dict_of_groups[groups_name] = cmd.get_object_list(groups_name)
    return dict_of_groups
cmd.extend("get_groups", get_groups)


def pic_every_group_alone(output_folder_path=op,bg="white",name="",width=0, height=0, dpi=300, ray=0, quiet=1,shots=[],group_map_name={}):
    if type(shots)== []:
        shots.append(cmd.get_view())
    cmd.bg_color(bg)
    all_obgects = cmd.get_names('objects', 0, '(all)')
    cmd.disable(" ".join(all_obgects))
    counter=0
    if type(shots) == type(8):
        shots_new=[]
        for i in range (shots):
            shots_new.append(cmd.get_view())
            cmd.turn("y", 360/shots)
        shots = shots_new
    groups = get_groups()
    for shot in shots:
        cmd.set_view(shot)
        for group in groups.keys():
            map_obj=group
            if group_map_name != {}:
                try:
                    map_obj = group_map_name[group[:6]]########the :6 is cutomize
                except:
                    continue
            cmd.enable(group)
            for obj in groups[group]:
                cmd.enable(obj)
            cmd.png(os.path.join(output_folder_path, "{}_{}_{}.png".format(counter,name,map_obj)), width=width, height=height, dpi=dpi, ray=ray, quiet=quiet)
            cmd.disable(group)
        counter+=1
cmd.extend("pic_every_group_alone", pic_every_group_alone)

#pic_every_group_alone(output_folder_path=op,bg="white",name="",width=0, height=0, dpi=500, ray=0, quiet=1,shots=8,group_map_name={})

# pic_every_obgect_alone(dpi=500,ray=0,shots=[
# (
#     -0.331894010,   -0.114593841,    0.936318994,
#     -0.942588925,    0.001509398,   -0.333930224,
#      0.036853522,   -0.993402839,   -0.108513713,
#      0.004345788,    0.004192799, -592.515319824,
#    269.653564453,  286.121704102,  289.544342041,
#   -68651.867187500, 69838.062500000,  -20.000000000 ),#side1
# (
#      0.466542512,   -0.131409124,    0.874668360,
#     -0.884422719,   -0.057306875,    0.463137716,
#     -0.010734294,   -0.989659965,   -0.142955169,
#      0.004345788,    0.004192799, -592.515319824,
#    269.653564453,  286.121704102,  289.544342041,
#   -68651.867187500, 69838.062500000,  -20.000000000 ),
# (
#      0.885688424,   -0.032991085,    0.463077843,
#     -0.464167148,   -0.082018077,    0.881932795,
#      0.008888676,   -0.996073544,   -0.087955020,
#     -0.003191456,    0.005846596, -592.619567871,
#    280.591522217,  279.909545898,  286.027557373,
#   -61727.367187500, 62913.562500000,  -20.000000000 ),
#
# (
#      0.196948603,    0.003109872,    0.980408072,
#     -0.976405323,   -0.089679852,    0.196428701,
#      0.088533387,   -0.995964348,   -0.014623457,
#      0.000797721,    0.000377506, -582.933593750,
#    270.122436523,  282.910339355,  285.624542236,
#   -66661.953125000, 67827.976562500,  -20.000000000 ),
# (
#      0.585947275,    0.055493567,    0.808445275,
#     -0.799701989,   -0.121566288,    0.587954402,
#      0.130908445,   -0.991028905,   -0.026850801,
#      0.001240727,    0.000220135, -582.955261230,
#    286.792297363,  259.866607666,  290.056854248,
#   -66661.953125000, 67827.976562500,  -20.000000000 ),
#
# (
#      0.641947806,   -0.764957428,   -0.052327219,
#     -0.766740978,   -0.640297055,   -0.046012402,
#      0.001694680,    0.069656983,   -0.997565091,
#      0.002011588,    0.002263963, -464.352569580,
#    288.595611572,  281.052062988,  290.548431396,
#   -66780.414062500, 67709.515625000,  -20.000000000 ),#close upper view
# (
#     -0.255235881,   -0.966741383,    0.015712244,
#     -0.955103695,    0.249566898,   -0.159625933,
#      0.150397196,   -0.055752721,   -0.987045467,
#      0.003665215,    0.002867922, -503.747619629,
#    270.070495605,  286.600952148,  290.841400146,
#   -68740.625000000, 69749.304687500,  -20.000000000 )#open upper view
#                               ])
#

#picture_roll(output_folder_path=op,bg="white",axis = "y",number_of_shots=8,angle=None,name="all_togther_no_singleton",width=0, height=0, dpi=500, ray=0, quiet=1)
#
# mapper={"cluster1": "6XEYKJ",
# "cluster2": "7JW0FG",
# "cluster3": "7KMHBA",
# "cluster4": "7LXYON",
# "cluster5": "7WOAGF",
# "cluster6": "7R40FD",
# "cluster7": "7PQZBA",
# "cluster8": "7BNVLH",
# "cluster9": "7DZXLH",
# "cluster10": "7AKDLH",
# "cluster11": "7RQ6KJ",
# "cluster12": "7FAELH",
# "cluster13": "7U2DLH",
# "cluster14": "8D6ZDC",
# "cluster15": "7RBVLH",
# "cluster16": "7S0ELH",
# "cluster17": "7RAQLH",
# "cluster18": "7SOBBC",
# "cluster19": "7UARGF",
# "cluster20": "7WK9ba"}
# mapper = {v:k for k,v in mapper.items()}
#
# pic_every_obgect_alone(output_folder_path=op,bg="white",name="all_togther",width=0, height=0, dpi=500, ray=0, quiet=1,object_map_name=mapper,
#  shots=8)
# op= "/run/user/29999/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/cov2_280822/regardent/analyzide_data/midi_bigger_cmpnd/pdb_abs_align/representive_1percluster_pics/s1"
# pic_every_obgect_alone(output_folder_path=op,bg="white",name="all_togther",width=0, height=0, dpi=500, ray=0, quiet=1,object_map_name=mapper,
#  shots=8)

#picture_roll(output_folder_path=op,bg="white",axis = "y",number_of_shots=8,angle=None,name="all_tother",width=0, height=0, dpi=500, ray=0, quiet=1)

# cmd.set_color( "cluster1", [0.8941176470588236, 0.10196078431372549, 0.10980392156862745])
# cmd.set_color(  "cluster2", [0.21568627450980393, 0.49411764705882355, 0.7215686274509804])
# cmd.set_color(  "cluster3", [0.30196078431372547, 0.6862745098039216, 0.2901960784313726])
# cmd.set_color(  "cluster4", [0.596078431372549, 0.3058823529411765, 0.6392156862745098])
# cmd.set_color(  "cluster5", [1.0, 0.4980392156862745, 0.0])
# cmd.set_color(  "cluster6", [1.0, 1.0, 0.2])
# cmd.set_color(  "cluster7", [0.6509803921568628, 0.33725490196078434, 0.1568627450980392])
# cmd.set_color(  "cluster8", [0.9686274509803922, 0.5058823529411764, 0.7490196078431373])
# cmd.set_color(  "cluster9", [0.10196078431372549, 0.8196078431372549, 1.0])
# cmd.set_color(  "cluster10", [0.3, 0.3, 0.3])
# cmd.set_color(  "cluster11",[0.4, 0.4, 0.4])
# cmd.set_color(  "cluster12", [0.8, 0.6, 1.0])
# cmd.set_color(  "cluster13",[0.7058823529411765, 0.6274509803921569, 0.39215686274509803])
# cmd.set_color(  "cluster14", [0.7019607843137254, 0.0, 0.5254901960784314])
# cmd.set_color(  "cluster15",[0.0, 0.2, 0.8])
# cmd.set_color(  "cluster16", [0.5,0.5,0.5])
# cmd.set_color(  "cluster17",[0.6,0.6,0.6])
# cmd.set_color(  "cluster18", [0.7,0.7,0.7])
# cmd.set_color(  "cluster19",[0.8,0.8,0.8])
# cmd.set_color(  "cluster20",[0.9,0.9,0.9])
#
# for obj in cmd.get_object_list('all'):
#     cmd.color(mapper[obj.split("_")[0]],obj+" and chain "+obj.split("_")[-1])
#     cmd.set_name(obj,mapper[obj.split("_")[0]]+"_"+obj.split("_")[0]+"_"+obj.split("_")[-1])
#
#


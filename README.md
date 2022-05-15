# color_peptides
The main purpose of this project is to color selected residues of your interested protein on a PDB with respect to a proper sequence alignment between your protein sequence and the sequence of the PDB file of your protein and generate a pictchers or gif, on pymol IDE.
## quick start
To do so you will need to change the last part of the python script of "coloring.py" show here. as it call the main function "color_pep_on_pdb()"
```
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
```
it will geneate 5 images as you can see in /example/outputs 
## little deeper

THIS IS THE MAIN FUNCTION OF THIS MODULE

it will color the list_of_residue_to_color on the pdb from your file path in pymol by the counting of the protein
seq that you gave by doing msa by muscle.
### parameters
name | type | explaination
--- | --- | ---
**pdb_file_path**            | (<ins>str</ins>)                   |   pdb file path can dowolond from any source of protein struchcer such as https://www.rcsb.org/  
**protein_seq**              | (<ins>str</ins>)                   | The sequence of the protein that you are working on                                              
**chains**                   |     (<ins>list</ins>)              | the chains on the pdb that your protein is based on i.e if it is a monomer so input the monomer chain id ["A"] if its a dimer that build by the chains "A","B" so input ["A","B"]. the order is relevant and should corespond to your protein_seq order !  
**list_of_residue_to_color** | (<ins>tupel of list or list</ins>) | tupel of ( list of residue to color base on the same order by your protein_seq )    or without the tupel for a single selection 
**color**                    |    (<ins>tupel of str or str</ins>)| tupel of ( str color as pymol names of color ) or without the tupel for a single selection 
**name_of_selection**        |   (<ins>tupel of str or str</ins>) |tupel of ( str name ) or without the tupel for a single selection. will be the name of pymol selectoin and the name of the output image. 
**output_folder_path**       |    (<ins>str</ins>)                | a path for all the output as image and video and muscle aligmnent. 
**muscle_folder_path**       |    (<ins>str</ins>)                | str the name of the folder that contain the binary file of muscle, the name of the binary file have to be muscle for linux or muscle.exe foe Windows, if muscle_folder_path = "" that mean that the muscle_folder_path is in the environmental path.
**together**                 | (<ins>bool</ins>)                  | spesify weather to take the image of all selection togher or for each selection make a seperate image. defult is False. for my opinion if their is not an overlap between the selection consider use together = True.
**len_movie**                |(<ins>int</ins>)                    | if you want to generate a mpg file(movie) that is the protein that roll against y axis so spedify a number of frame you would like to have 240 is 8 second len of movie
**dpi**                      | (<ins>int</ins>)                   | the dpi of the image 150 is fine for a scatch consider dpi=300 for an article view
**ray**                      |(<ins>int 0 or 1 </ins>)            | is to make the image ray, it add a value to the 3D texture it is nice for an article but it is a time  expensive

import random

prot_list = [['M1','C']]
locations_list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
aa_list = ['VAL', 'ILE', 'LEU', 'GLU', 'GLN', 'ASP', 'ASN', 'HIS', 'TRP', 'PHE', 'TYR',
'ARG', 'LYS', 'SER', 'THR', 'MET', 'ALA', 'GLY', 'PRO', 'CYS']

locations = random.choices(locations_list,weights=None,k=9) 
aa = random.choices(aa_list,weights=None,k=9)

for prot in prot_list:
    #Initialize
    # load EQ MART1 (3qdj)
    cmd.select(prot[0])

    # mutagenesis mode
    cmd.wizard("mutagenesis")
    cmd.do("refresh_wizard")
    
    # create 9 objects 
    for obj in range(1,10):
        cmd.create('obj_%s' % obj, prot[0])
    # Mutate
    cmd.get_wizard().set_mode(aa[0])
    
    #for obj_1 (random aa mutation(s) in random location(s))   
    cmd.get_wizard().do_select("/obj_1//%s/%d" % (prot[1], locations[0]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[0], aa[0]),'obj_1')
    
    #for obj_2 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[1])
    cmd.get_wizard().do_select("/obj_2//%s/%d" % (prot[1], locations[1]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[1], aa[1]),'obj_2')
   
    #for obj_3 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[2])
    cmd.get_wizard().do_select("/obj_3//%s/%d" % (prot[1], locations[2]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[2], aa[2]),'obj_3')
   
    #for obj_4 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[3])
    cmd.get_wizard().do_select("/obj_4//%s/%d" % (prot[1], locations[3]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[3], aa[3]),'obj_4')

    #for obj_5 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[4])
    cmd.get_wizard().do_select("/obj_5//%s/%d" % (prot[1], locations[4]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[4], aa[4]),'obj_5')

    #for obj_6 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[5])
    cmd.get_wizard().do_select("/obj_6//%s/%d" % (prot[1], locations[5]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[5], aa[5]),'obj_6')

    #for obj_7 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[6])
    cmd.get_wizard().do_select("/obj_7//%s/%d" % (prot[1], locations[6]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[6], aa[6]),'obj_7')

    #for obj_8 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[7])
    cmd.get_wizard().do_select("/obj_8//%s/%d" % (prot[1], locations[7]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[7], aa[7]),'obj_8')

    #for obj_9 (random aa mutation(s) in random location(s))
    cmd.get_wizard().set_mode(aa[8])
    cmd.get_wizard().do_select("/obj_9//%s/%d" % (prot[1], locations[8]))
    cmd.frame(1)
    cmd.get_wizard().apply()
    cmd.save("/Users/zrollins/Documents/Documents/DMF5_MART1/%s_%s_%s_%s.pdb" % (prot[0], prot[1], locations[8], aa[8]),'obj_9')

 # Done
    cmd.set_wizard()
    if len(prot_list) == 1:
        continue
    else:
        cmd.reinitialize()  
   

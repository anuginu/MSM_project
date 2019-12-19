from villin_features import *
from ga_functions import *
from calculate_gmrq_score import *
dih=Get_dihedral_features_villin()
#cont=Get_contacts_features_villin()
#comb=Get_combined_features_villin()

score_spec_dih=spec_score(dih)
score_lap_dih=Laplacian_score(dih)
score_mcfs_dih=mcfs_score(dih)
score_ndfs_dih=ndfs_score(dih)
score_udfs_dih=udfs_score(dih)
score_spec_feature=[]
for i in range(10,len(score_spec_dih[1]),10):
  dih_gmrq=calculate_scor(score_spec_dih[1][0:i],dih)
  score_spec_feature.append(dih_gmrq)

print(score_spec_feature) 
f=open("villin_dih_spec_ranking.txt","a+")  
f.write(str(score_spec_feature))
f.close()
score_lap_feature=[]
for i in range(10,len(score_lap_dih[1]),10):
  dih_lap=calculate_scor(score_lap_dih[1][0:i],dih)
  score_lap_feature.append(dih_lap)
f=open("villin_dih_lap_ranking.txt","a+")
f.write(str(score_lap_feature))
f.close()



score_mcfs_feature=[]
for i in range(10,len(score_mcfs_dih[1]),10):
  dih_mcfs=calculate_scor(score_mcfs_dih[1][0:i],dih)

  score_mcfs_feature.append(dih_mcfs)
f=open("villin_dih_mcfs_ranking.txt","a+")
f.write(str(score_mcfs_feature))
f.close()


score_ndfs_feature=[]
for i in range(10,len(score_ndfs_dih[1]),10):
  dih_ndfs=calculate_scor(score_ndfs_dih[1][0:i],dih)

  score_ndfs_feature.append(dih_ndfs)
f=open("villin_dih_ndfs_ranking.txt","a+")
f.write(str(score_ndfs_feature))
f.close()



score_udfs_feature=[]
for i in range(10,len(score_udfs_dih[1]),10):
  dih_udfs=calculate_scor(score_udfs_dih[1][0:i],dih)

  score_udfs_feature.append(dih_udfs)
f=open("villin_dih_udfs_ranking.txt","a+")
f.write(str(score_udfs_feature))
f.close()


from villin_features import *
from ga_functions import *
from calculate_gmrq_score import *
cont=Get_combined_features_villin()

score_spec_cont=spec_score(cont)
score_spec_feature=[]
for i in range(10,len(score_spec_cont[1]),100):
  cont_gmrq=calculate_scor(score_spec_cont[1][0:i],cont)
  score_spec_feature.append(cont_gmrq)

print(score_spec_feature) 
f=open("villin_cont_spec_ranking.txt","a+")  
f.write(str(score_spec_feature))
f.close()  
score_lap_cont=Laplacian_score(cont)

score_lap_feature=[]
for i in range(10,len(score_lap_cont[1]),100):
  cont_gmrq=calculate_scor(score_lap_cont[1][0:i],cont)
  score_lap_feature.append(cont_gmrq)

print(score_lap_feature) 
f=open("villin_cont_lap_ranking.txt","a+")  
f.write(str(score_lap_feature))
f.close()  



score_mcfs_cont=mcfs_score(cont)

score_mcfs_feature=[]
for i in range(10,len(score_mcfs_cont[1]),100):
  cont_gmrq=calculate_scor(score_mcfs_cont[1][0:i],cont)
  score_mcfs_feature.append(cont_gmrq)

print(score_mcfs_feature) 
f=open("villin_cont_mcfs_ranking.txt","a+")  
f.write(str(score_mcfs_feature))
f.close()  

score_ndfs_cont=ndfs_score(cont)

score_ndfs_feature=[]
for i in range(10,len(score_ndfs_cont[1]),100):
  cont_gmrq=calculate_scor(score_ndfs_cont[1][0:i],cont)
  score_ndfs_feature.append(cont_gmrq)

print(score_ndfs_feature) 
f=open("villin_cont_ndfs_ranking.txt","a+")  
f.write(str(score_ndfs_feature))
f.close()  







score_udfs_cont=udfs_score(cont)


score_udfs_feature=[]
for i in range(10,len(score_udfs_cont[1]),100):
  cont_gmrq=calculate_scor(score_udfs_cont[1][0:i],cont)
  score_udfs_feature.append(cont_gmrq)

print(score_udfs_feature) 
f=open("villin_cont_udfs_ranking.txt","a+")  
f.write(str(score_udfs_feature))
f.close()  

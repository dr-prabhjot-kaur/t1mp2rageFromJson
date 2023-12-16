# t1mp2rageFromJson

1. 
Either pull from dockerhub


docker pull prabhjotkaur1991/t1mp2rage:latest 

OR

Build it from github repo


git clone https://github.com/dr-prabhjot-kaur/t1mp2rageFromJson.git


docker build username/t1mp2rage:latest

2. Run the docker file to estimate t1 map from mp2rage images.



requirement: UNI nii image, INV1 json, INV2json



docker run -v /home:/home prabhjotkaur1991/t1mp2rage:latest sh run_EstimateT1fromMP2RAGE1.sh /opt/mcr/R2023a  '/home/ch239656/work/t1w_MP2RAGE_UNI/t1w_MP2RAGE_UNI.nii.gz' '/home/ch239656/work/t1w_MP2RAGE_INV1/t1w_MP2RAGE_INV1.json' '/home/ch239656/work/t1w_MP2RAGE_INV2/t1w_MP2RAGE_INV2.json' 0.0064 30 150 1 /home/output.nii.gz

#help:


#TRFLASH='0.0064'


#NZslices1='30'


#NZslices2='150'


#eff='1'

#outputfilename='/home/output.nii.gz'

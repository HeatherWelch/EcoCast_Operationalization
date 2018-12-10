# Code library to make operational EcoCast - a dynamic ocean management tool for fisheries bycatch
### Code authors: Heather Welch (UCSC, NOAA), Elliott Hazen (NOAA), Kylie Scales (University of the Sunshine Coast)

<img src="inst/imgs/logo.png?raw=True" width="400">

### Relevant manuscripts: 
**Welch et al. 2018.** "Practical considerations for operationalizing dynamic management tools." Journal of Applied Ecology 54: 1415-1428. doi: 10.1111/1365-2664.13281.  
**Hazen et al. 2018.** “A dynamic ocean management tool to reduce bycatch and support sustainable fisheries.” Science Advances 4: eaar3001.    
**Hazen et al. 2017.** "WhaleWatch: a dynamic management tool for predicting blue whale density in the California Current." Journal of Applied Ecology 54: 1415-1428.  

### Description of scripts
**1. Get_Env_Data_A.R** (run once at beginning of the day): Get data sequence number one : Create final and temporary directories, acquire all static variables.  
**2. Get_Env_Data_B.R** (run multiple times during the day): Get data sequence number two : See which dynamic variables are missing. If none are missing, run EcoCast. If variables are missing, attempt to download missing variables. See which dynamic variables are still missing after download attempt. If none are missing, run EcoCast.  
**3. Get_Env_Data_C.R** (run once at end of the day): Get data sequence number three : Evaluates the most recent version of environmental data available, and then runs EcoCast. Script will not overwrite pre-existing final products (i.e. if final products were created by Get_Env_Data_B.R).  


For EcoCast development, see Elliott Hazen's github repository: https://github.com/elhazen/EcoCast-SciAdv

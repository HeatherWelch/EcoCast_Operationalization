# Code library to make operational EcoCast - a dynamic ocean management tool for fisheries bycatch
### Code authors: Heather Welch (UCSC, NOAA), Elliott Hazen (NOAA), Kylie Scales (University of the Sunshine Coast)

<img src="inst/imgs/logo.png?raw=True" width="400">

### Relevant manuscripts: 
Hazen et al. “An Eco-Informatic solution to ocean bycatch.” In review.  
Welch et al. "Practical considerations for operationalizing dynamic management tools." In prep.  
Hazen et al. 2017. "WhaleWatch: a dynamic management tool for predicting blue whale density in the California Current." Journal of Applied Ecology 54: 1415-1428.

### Description of scripts
1. Get_Env_Data_A.R (run once at beginning of the day): Get data sequence number one : Create final and temporary directories, acquire all static variables.
2. Get_Env_Data_B.R (run multiple times during the day): Get data sequence number two : See which dynamic variables are missing. If none are missing, run EcoCast. If variables are missing, attempt to download missing variables. See which dynamic variables are still missing after download attempt. If none are missing, run EcoCast.
3. Get_Env_Data_C.R (run once at end of the day): Get data sequence number three : Evaluates the most recent version of environmental data available, and then runs EcoCast. Script will not overwrite pre-existing final products (i.e. if final products were created by Get_Env_Data_B.R).

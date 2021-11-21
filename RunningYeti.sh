ssh tcline@yeti.cr.usgs.gov

#Enter active directory password
#Account norock
#User tcline


#Move data
scp test.txt tcline@yeti.cr.usgs.gov:~/

mkdir Annual_NorWest_Data
mv test.txt data/

scp -r data/ username@yeti.cr.usgs.gov:~/

#Navigate to where you would like to transfer files from
#in this case Cline_USGS hard drive
scp -r Annual_NorWest_Data/Regions/Spokoot.ssn tcline@yeti.cr.usgs.gov:~/Annual_NorWest_Data/Regions
scp Spokoot.DayMetSSN.RDS tcline@yeti.cr.usgs.gov:~/Annual_NorWest_Data

#module avail
module list
module load R/4.1.1
module load gdal/3.1.0
module load proj/7.0.1

##Install rgdal
##first load gdal and proj modules
##make sure 'sp' is installed in R
#wget https://cran.r-project.org/src/contrib/rgdal_1.5-27.tar.gz
#R CMD INSTALL rgdal_1.5-27.tar.gz --configure-args="--with-proj-include=/opt/ohpc/pub/usgs/libs/gnu8/proj/7.0.1/include --with-proj-lib=/opt/ohpc/pub/usgs/libs/gnu8/proj/7.0.1/lib --with-proj-share=/opt/ohpc/pub/usgs/libs/gnu8/proj/7.0.1/share/proj"

#module load jags/4.3.0

R #open R

#install.packages('lme4')
#install.packages('sp')
#install.packages('SSN')
#install.packages('here')
#install.packages('foreign')
#install.packages('dplyr')
#install.packages('tidyr')
#install.packages('doParallel')
#install.packages('foreach')

q() #quit R

#module unload R/4.0.1

#git clone https://code.chs.usgs.gov/sas/arc/Getting_Started_Yeti_Examples.git
#ls
#cd ~
#rm -r Annual_NorWest
#git clone https://github.com/timothycline/Annual_NorWest.git

git pull

#git reset --hard https://github.com/timothycline/Annual_NorWest.git/master

cd Annual_NorWest
#cat foreach.slurm 
#vi foreach.slurm

#sbatch test.slurm #Submit Job
#cat SpokootDayMet.out

#sbatch DayMetRefit_Parallel.slurm #Submit Job #5840483
#cat DayMetRefitParallel.out
sbatch Spokoot.slurm
cat SpokootSSN_Refit.out

sbatch Clearwater.slurm
cat ClearwaterSSN_Refit.out

sbatch SnakeBear.slurm
sbatch MissouriHW.slurm
sbatch Spokoot.slurm
sbatch UpMissMarias.slurm
sbatch 

sidle #Check Status


scancel <jobid>
scancel 5840479
scancel -u <tcline>


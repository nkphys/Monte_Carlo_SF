###### RUN With: bash BBB.sh ######################

cwd=$PWD
datafile=input.inp
#sub=job.sge
exe=sf


  for Beta in 50 #80 100 120 150 #150 290 #{40..400..40}   
   do
    BetaDIR=$cwd'/Beta-'$Beta
     mkdir $BetaDIR
     cp $exe $BetaDIR
     cd $BetaDIR
     echo $BetaDIR

#Make the datafile
        rm $datafile
#        rm $sub

	echo "Xsite=8" >> $datafile
	echo "Ysite=8" >> $datafile
	echo "Orbitals=3" >> $datafile
	echo "Temperature=$Beta" >> $datafile
	echo "Fill=0.666666667" >> $datafile
	echo "MaxMCsteps=1000" >> $datafile
	echo "RandomSeed=7915" >> $datafile
	echo "J_NN=0.012" >> $datafile
	echo "J_NNN=0.008" >> $datafile
  
#make the job script
#	echo "#$ -N Beta"$Beta"CHEM"$CHEM >> $sub
#	echo "#$ -q medium_sigma*" >> $sub
#	echo "#$ -cwd" >> $sub
#	echo "date" >> $sub
#	echo "./$exe " >> $sub
#	echo "date" >> $sub
#	#qsub $sub    #to run on laptop comment this line
      time ./$exe $datafile > out.txt &    ## to run on laptop uncomment

done 
cd $cwd

#!/bin/sh


# download CRPropa3-data repository
# this will provide additional tools to compute the interaction rate tables
if [ -d "CRPropa3-data" ] 
then
	echo "CRPropa3-data already exists. Continuing..." 
else
	echo "CRPropa3-data already exists. Downloading it..."
	git clone https://github.com/CRPropa/CRPropa3-data.git
fi

# directory where data will be stored
dataDir="../build/data"

# folders with tabulated interactions
folders=("PairProduction", "InverseComptonScattering")


# compute the tables
echo "Generating the interaction tables..."
$PYTHON_EXECUTABLE computeElectromagneticInteractions.py

# copy to the installation directory (assumed to be build/data)
if [ -d ${dataDir} ]
then

	interactionDataDir="${dataDir}/${folder}"
	if ! [ -d ${interactionDataDir} ]
	then
		mkdir ${interactionDataDir}
	fi

	for folder in ${folders[@]}
	do

		if [ -d "../build/data/${folder}" ]
		then
			cp -r ../data/${folder}/* ../build/data/${folder}/.
			cp -r ../data/${folder}/* ../build/share/alpinist/${folder}/.
		fi

	done

	echo "Files copied to the appropriate folder." 
else
	echo "Installation path is not standard." 
	echo "You should mannually copy the data files from ../data/ to the folder where ALPinist was installed."
fi
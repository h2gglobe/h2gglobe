##
## Common options
##

[[ -z $version  ]] && [[ -f version.sh ]] && . version.sh

base_storedir="./datastore"
## storeremote="/castor/cern.ch/user/c/cmshgg/reduced"
storeremote="/store/group/phys_higgs/cmshgg/reduced/${group}"

export storeremote

[[ -f $(whoami)_setup.sh ]] && . $(whoami)_setup.sh

if [ -n "${BATCH+x}" ]; then
    [[ -d ${base_storedir} ]] && rm -rf ${base_storedir}
    mkdir ${base_storedir}
    chmod -R 775 ${base_storedir}
else 
    [[ ! -d ${base_storedir} ]] && mkdir ${base_storedir}
    if [[ -n ${storeremote} ]]; then
	if ( ! echo ${storeremote} | egrep '(castor|/store)' > /dev/null ) && 
	    ( ! mount | grep ${storeremote} | grep ${USER} > /dev/null ); then 
	    sshfs ${storeremote} ${storedir} -o nonempty
	fi
    fi
fi

storedir=${base_storedir}/${version}

[[ ! -d ${storedir} ]] && mkdir $storedir

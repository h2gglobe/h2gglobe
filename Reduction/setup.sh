##
## Common options
##

[[ -z $version  ]] && [[ -f version.sh ]] && . version.sh

base_storedir="./datastore"
storedir="${base_storedir}"
storeremote="/castor/cern.ch/user/c/cmshgg/reduced"

[[ -f $(whoami)_setup.sh ]] && . $(whoami)_setup.sh

remoteneeded="inputs"

if [ -n "${BATCH+x}" ]; then
	[[ -d ${base_storedir} ]] && rm -rf ${base_storedir}
	mkdir ${base_storedir}
	chmod -R 775 ${base_storedir}
	if [[ -n ${storeremote} ]] && ( ! echo ${storeremote} | grep castor > /dev/null ); then
  	    for d in $remoteneeded; do
		rsync -avz ${storeremote}/${d} ${base_storedir}
	    done	
	fi
else 
    [[ ! -d ${storedir} ]] && mkdir ${storedir}
    if [[ -n ${storeremote} ]]; then
	if ( ! echo ${storeremote} | grep castor > /dev/null ) && 
	    ( ! mount | grep ${storeremote} | grep ${USER} > /dev/null ); then 
	    sshfs ${storeremote} ${storedir} -o nonempty
	fi
    fi
fi

storedir=${storedir}/${version}

[[ ! -d ${storedir} ]] && mkdir $storedir

hset () {
  eval hash"$1$2"='$3'
}

hget () {
  eval echo '${hash'"$1$2"'#hash}'
}

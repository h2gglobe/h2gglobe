if [ -n "${BATCH+x}" ]; then
    if ( ! echo ${storeremote} | grep castor  ); then
	rsync -avz ${storedir}/ ${storeremote}/${version}
    else 
	rfmkdir ${storeremote}/${version}
	cd ${storedir}
	for d in $(find -type d); do
	    rfmkdir ${storeremote}/${version}/$d;
	    for f in $(find $d -maxdepth 1 -type f); do
		rfcp $f ${storeremote}/${version}/$d &
	    done
	done
    fi
fi


wait


sha1=""
errs=""
if [ -n "${BATCH+x}" ]; then
    ls -R ${storedir}

    if ( ! echo ${storeremote} | egrep '(castor|store)'  ); then
	rsync -avz ${storedir}/ ${storeremote}/${version}
    else 
	if( echo ${storeremote} | grep 'castor' ); then
	    mkdir=rfmkdir
	    cp=rfcp
	else
	    mkdir=cmsMkdir
	    cp=cmsStage
	fi
	$mkdir ${storeremote}/${version}
	cd ${storedir}
	for d in $(find -type d); do
	    $mkdir ${storeremote}/${version}/$d;
	    for f in $(find $d -maxdepth 1 -type f); do
		sha1="$sha1\n$(sha1sum $f)"
		$cp -f $f ${storeremote}/${version}/$d
		exstat=$?
		if [[ "$exstat" != 0 ]]; then
		    errs="$errs $f($exstat)"
		fi
		
	    done
	done
    fi
fi


wait

echo -e $sha1

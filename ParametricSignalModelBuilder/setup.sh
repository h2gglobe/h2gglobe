#!/bin/sh

# run this as 
#
#   source setup.sh
#

if [[ $(hostname) == uaf* ]] ; then
  CONFIG_HOSTNAME=uaf
elif [[ $(hostname) == lxplus* ]] ; then
  CONFIG_HOSTNAME=lxplus
else
  CONFIG_HOSTNAME=$(hostname -f)
fi

if [ -e setup-${CONFIG_HOSTNAME}.sh ]; then
  echo "using setup-${CONFIG_HOSTNAME}.sh"
  source setup-${CONFIG_HOSTNAME}.sh
else
  echo "don't know how to set up the environment on this host (specific configuration file not found)"
fi
#!/bin/sh

SSLVERSION='0.9.8'
INSTALL_PATH='/usr/local'

ln_ssl_crypto()
{
    if [ -d $1 ]; then
        cd $1
        echo changed directory to $1
        if [ ! -f libssl.so.4 ]; then
            ln -s libssl.so.$SSLVERSION libssl.so.4
        fi
        if [ ! -f libcrypto.so.4 ]; then
            ln -s libcrypto.so.$SSLVERSION libcrypto.so.4
        fi
    fi
}

if [ "$?" -gt "0" ]
then
echo the installation ended
echo error: no archive given as argument or file not found
echo e.g. "sh install.sh NRGsuite.tar.gz"
exit 1
fi

echo cleaning...
rm -rf ${INSTALL_PATH}/NRGsuite

echo un-archiving...
tar -zxvf $1

echo copying files...
cp -a NRGsuite ${INSTALL_PATH}

if [ "$?" -gt "0" ]
then
echo the installation ended
echo error: could not copy NRGsuite folder
exit 1
fi

echo changing permissions...
chmod a+x ${INSTALL_PATH}/NRGsuite/FlexAID/WRK/FlexAID
chmod a+x ${INSTALL_PATH}/NRGsuite/FlexAID/WRK/Process_Ligand
chmod a+x ${INSTALL_PATH}/NRGsuite/GetCleft/WRK/GetCleft
chmod a+x ${INSTALL_PATH}/NRGsuite/GetCleft/WRK/volume_calc

echo creating symbolic links...
cd ${INSTALL_PATH}/NRGsuite/FlexAID/WRK/libs
echo changed directory to ${INSTALL_PATH}/NRGsuite/FlexAID/WRK/libs
ln -s libinchi.so.0.4.1 libinchi.so.0
ln -s libinchi.so.0 libinchi.so
ln -s libopenbabel.so.4.0.2 libopenbabel.so.4
ln -s libopenbabel.so.4 libopenbabel.so

ln_ssl_crypto /lib/
ln_ssl_crypto /lib32/
ln_ssl_crypto /usr/lib/
ln_ssl_crypto /usr/lib32/

echo installation successful
echo done.

echo "" >> "$HOME/.pymolrc"
echo "import os" >> "$HOME/.pymolrc"
echo "os.environ['NRGSUITE_INSTALLATION'] = \"${INSTALL_PATH}/NRGsuite\"" >> "$HOME/.pymolrc"

exit 0

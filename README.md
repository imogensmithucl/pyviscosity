# pyviscosity


# SSL certificates for MacOS

In some cases, to use the urlib library within cirpy, one first need to update the SSL certificates.
There are scripts within the /Library directory, find them with:
    mdfind -name "Install Certificates"
    /Applications/Python 3.12/Install Certificates.command
    /Applications/Python 3.11/Install Certificates.command
    /Applications/Python 3.7/Install Certificates.command
    /Applications/Python 3.8/Install Certificates.command


Then execute those scripts (some might be failing, if the corresponding python version is not installed, but this is not a problem)

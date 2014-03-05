INSTALLATION (WINDOWS):

MARMITES is programmed in Python (www.python.org).

STEP1
Install first Python version 2.7.2 (download at http://www.python.org/download/releases/2.7.2/).
You also need to install the following packages:
- pylab (i.e. numpy+scipy) (www.scipy.org)
- matplotlib               (matplotlib.sourceforge.net)
- wxPython                 (www.wxpython.org)
- h5py                     (code.google.com/p/h5py/)
- pywin32                  (sourceforge.net/projects/pywin32/)
Another alternative is to install Python(x,y), see http://code.google.com/p/pythonxy/ and download the installation file at http://ftp.ntua.gr/pub/devel/pythonxy/Python(x,y)-2.7.2.1.exe.
Note that you may have to update some packages to the version used by MARMITES.
It is also recommended to install a development environment such as PyScripter (http://code.google.com/p/pyscripter/) or Spyder (this last one is already integrated in Python(x,y)).
As the MARMITES and MODFLOW data are stored in HDF5 format, you may also want to install HDFView (http://www.hdfgroup.org/hdf-java-html/hdfview/) to visualize and export these data.

STEP2
Install TortoiseHg (download at http://tortoisehg.bitbucket.org/)
Install TortoiseSVN (download at http://tortoisesvn.net/)

STEP3
- create a MARMITES folder in your local Python (typically C:/Python27) in the folder %pythonpath%\Lib\site-packages
- right click on the MARMITES folder, select TortoiseHg and Clone
- in the Source box, type https://frances.alain@code.google.com/p/marmites/ and press Clone
- copy the file MARMITES.pth (currently under MARMITES folder) into the folder %pythonpath%\Lib\site-packages

STEP4
- create a FloPy folder in your local Python (typically C:/Python27) in the folder %pythonpath%\Lib\site-packages
- right click on the FloPy folder, select SVN Checkout...
- in the URL of repository box, type https://flopy.googlecode.com/svn and press Ok
- copy the file FloPy.pth (currently under MARMITES folder) into the folder %pythonpath%\Lib\site-packages

STEP5
- install MODFLOW-NWT, download at http://water.usgs.gov/nrp/gwsoftware/modflow_nwt/MODFLOW-NWT_1.0.8.zip
- unzip in any folder (typically C:/MODFLOW-NWT_1.0.8)
To link the soil reservoir to the groundwater reservoir, the UZF1 package is used (see http://code.google.com/p/marmites/source/browse/doc/MARMITES_schema_201110.png).
Documentation of the UZF1 package can be found at http://pubs.usgs.gov/tm/2006/tm6a19.

STEP6
- prepare the required input as in the example located in the doc folder (if not updated contact me at frances.alain@gmail.com).
The input are organized as follow: one head folder/workspace MM_ws with the MMunsat input and two sub folders, MF_ws with MODFLOW input and MMsurf_ws with MARMITES surface input. Each one of the 3 folders contains an .ini file that defines all the parameters of the MARMITES and MODFLOW modules.
Generally you will have to update: (i) the path of the MMsurf_ws and MF_ws in __inputMF.ini, as well as the name of the files with the data; (ii) the MODFLOW path and name in the __inputMF.ini.
See the example under MARMITES/doc/example. Open the __input*.ini files, they contain instructions and explanation about the meaning of the parameters. For MODFLOW the name of the variables is kept as in the MODFLOW manual.
- edit the file startMM_fn.txt (on desktop) and update it with the correct path and name of your main MM ini file.
- start MARMITES with startMARMITES_v3.py

A manual is under construction, as well as several applications and papers.
Images output can be found in MARMITES/doc/imgs. You can run a simple synthetic model using the data stored in MARMITES/doc/example.

REFERENCES:
Francés, A.P., Berhe, E., Lubczynski, M.W., 2010. Spatio - temporal groundwater recharge assessment using a lumped - parameter distributed model of the unsaturated zone, pyEARTH-2D, In: Geophysical Research Abstracts, 12(2010) EGU2010-6627-2, EGU general assembly 2010. 2 p.

Francés, A.P., Reyes-Acosta, J.L., Balugani, E., van der Tol, C., Lubczynski, M.W., 2011. Assessment of catchment water balance using distributed and transient coupled models of the unsaturated and saturated zones. Presented at ModelCare 2011 : models : repositories of knowledge, 18-22 September 2011, Leibzig, Germany. 2 p.

Francés, A.P., Reyes-Acosta, J.L., Balugani, E., van der Tol, C., Lubczynski, M.W., 2011. Towards an improved assessment of the water balance at the catchment scale : a coupled model approach. In: Estudios en la zona no saturada del suelo : volumen X : ZNS11 proceedings, 19-21 October 2011, Salamanca, Spain : e-book / editor J.M. Fernández, , N.S. Martin. - Salamanca : Universidad de Salamanca, 2011. - 370 p. ISBN 978-84-694-6642-1. pp. 321-326.
Available at: http://www.zonanosaturada.com/zns11/publications/p321.pdf


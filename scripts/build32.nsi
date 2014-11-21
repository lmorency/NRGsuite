# define the name of the installer
Outfile "NRGsuite_Win32.exe"

# define the directory to install to, the desktop in this case as specified  
InstallDir "$PROGRAMFILES\NRGsuite"

PageEx directory
  DirText "" "" "" ""
PageExEnd

Section install
  # define the output path for this file
  SetOutPath $INSTDIR

  # define what to install and place it in the output path
  File /r /x *.py~ /x *.pyc /x scripts "C:\Program Files\NRGsuite\*"
SectionEnd

PageEx instfiles
PageExEnd

Section writerc
  FileOpen $0 "$PROFILE\pymolrc.pml" a
  FileWrite $0 `import os$\r$\n`
  FileWrite $0 `os.environ['NRGSUITE_INSTALLATION'] = "$INSTDIR"$\r$\n`
  FileClose $0
SectionEnd
